from argschema import ArgSchemaParser
import os
import shutil
import sys
import subprocess
import time
import fnmatch
from pathlib import Path

import numpy as np

from ...common.utils import catGT_ex_params_from_str
from ecephys_spike_sorting.scripts.helpers import SpikeGLX_utils

def call_TPrime(args):

    # Run TPrime on a "standard" multiprobe + NI, using NP 1.0 or 2.0, with run
    # folder and probe folders
    # inputs:
    # full path to Tprime executable
    # parameters for "to stream" extracted edges ("to stream" = the reference stream for the data set)
    # parameters for NI sync edges
    
    # bNPY, if True, no text files of spike times will be written
    bNPY = True;
    if bNPY:
        outSuffix = '.npy'
    else:
        outSuffix = '.txt'

    print('ecephys spike sorting: TPrime helper module')
    start = time.time()
    
    # build paths to the input data for TPrime
    first_gate, last_gate = SpikeGLX_utils.ParseGateStr(args['catGT_helper_params']['gate_string'])
    catGT_dest = args['directories']['extracted_data_directory']
    run_name = args['catGT_helper_params']['run_name'] + '_g' + str(first_gate)
    run_dir_name =  args['tPrime_helper_params']['catGT_out_tag'] + '_' + run_name
    prb_dir_prefix = run_name + '_imec'
    
    run_directory = os.path.join( catGT_dest, run_dir_name ) # extracted edge files for aux data reside in run directory
    sync_period = args['tPrime_helper_params']['sync_period']
    
    sort_out_tag = args['tPrime_helper_params']['sort_out_tag']
    output_pat = 'imec*_' + sort_out_tag
    sort_out_list = list()
    
    # check for presence of an fyi file, indicating run with catgt 3.0 or later
    fyi_path = run_directory.replace('\\', '/') + '/' + run_name + '_all_fyi.txt'
    all_fyi_exists = Path(fyi_path).is_file()

    if args['tPrime_helper_params']['catGT_out_tag'] == 'supercat':
        # only use offsets for runs that have been processed with supercat
        # output from catgt should not have discontinuites.
        # note that to use offsets for non-supercat data, the offsets would
        # need to be compiled in catgt_helper (just like the fyi_all)
        offset_path =  run_directory.replace('\\', '/') + '/' + run_name + '_sc_offsets.txt'        
        offset_exists = Path(offset_path).is_file()
    else:
        offset_exists = False
    
    if not all_fyi_exists:
        # check for an fyi file -- assume that CatGT was run directly from
        # a batch file for all probes, and use that, but also print message
        print('No _all_fyi.txt file found.')
        print('Checking for _fyi.txt file from CatGT run outside ecephys pipeline.')
        fyi_path = run_directory.replace('\\', '/') + '/' + run_name + '_fyi.txt'
        fyi_exists = Path(fyi_path).is_file()

    if all_fyi_exists or fyi_exists:
        
        toStream_params = args['tPrime_helper_params']['toStream_sync_params']
        toStream_js, toStream_ip = parse_stream(toStream_params)
        toStream_id = (toStream_js, toStream_ip)
        toStream_path, from_list, from_list_ids, events_list, from_stream_index, out_list, all_list \
            = parse_catgt_fyi(fyi_path, toStream_id)
            
        
        if toStream_js == 2:
            # if toStream is an imec probe, create the file of spike times in sec for all output files
            prb_dir = prb_dir_prefix + str(toStream_ip)
            prb_path =  Path(os.path.join(run_directory, prb_dir)) 
            # get list of output directories -- can be multiples when splitting shanks
            ks_outdir_iter = prb_path.glob(output_pat)
            for ks_outdir in ks_outdir_iter: 
                sort_out_list.append(ks_outdir)
                st_file = os.path.join(run_directory, prb_dir, ks_outdir, 'spike_times.npy')
                # convert to seconds; if bNPY = True, returned file is an npy file
                # otherwise, text.
                toStream_events_sec = spike_times_npy_to_sec(st_file, 0, bNPY)
                # create a copy with name = spike_times_sec.adj
                shutil.copyfile(toStream_events_sec, os.path.join(run_directory, prb_dir, ks_outdir, 'spike_times_sec_adj.npy'))
                # if data was saved as text, also save as npy
                if not bNPY:
                    spike_times_sec_to_npy(toStream_events_sec)
                        
        # loop over the from_list_ids; for any that are probes, need to create 
        # files of spike times, and append to event_list, from_stream_index, and out_list
        for i, id in enumerate(from_list_ids):
            if id[0] == 2:      # imec stream                
                prb_dir = prb_dir_prefix + str(id[1])
                prb_path = Path(os.path.join(run_directory, prb_dir))                
                # get list of output directories -- can be multiples when splitting shanks
                ks_outdir_iter = prb_path.glob(output_pat)
                for ks_outdir in ks_outdir_iter:
                    sort_out_list.append(ks_outdir)
                    st_file = os.path.join(run_directory, prb_dir, ks_outdir, 'spike_times.npy')
                    # convert to seconds; if bNPY = True, returned file is an npy file
                    # otherwise, text.
                    st_file_sec = spike_times_npy_to_sec(st_file, 0, bNPY)
                    from_stream_index.append(i)
                    events_list.append(st_file_sec)                
                    # build path for output spike times text file
                    out_name = 'spike_times_sec_adj' + outSuffix
                    out_file = os.path.join(run_directory, prb_dir,ks_outdir, out_name)
                    out_list.append(out_file)
                         
    else:
        # Must be running with data from older CatGT with no fyi file
        
        print('No fyi file found.')
        print('You can generate a new fyi file by running an extract only pass on the CatGT output.')
        print('Typical extract only command line also see CatGT Readme under: -t=cat defer extraction to a later pass ')
        print('-dir=<catGT_dest> -run=catgt_<run_name> -g=0 -t=cat -ap -ni -prb=0 -prb_fld -no_tshift ^ ' )
        print('-ap -ni -prb=0 -prb_fld -no_tshift ^')
        print('<ni extraction string for CatGT 4.0> ^')
        print('dest=<catGT_dest>')    


    print('toStream:')
    print(toStream_path)
    print('fromStream')
    for fp in from_list:
        print(fp)
    print('event files')
    for i, ep in enumerate(events_list):
        print('index: ' + repr(from_stream_index[i]) + ',' + ep)
    print('output files')
    for op in out_list:
        print(op)
        
    # path to the 'runit.bat' executable that calls TPrime.
    # Essential in linux where TPrime executable is only callable through runit
    if sys.platform.startswith('win'):
        exe_path = os.path.join(args['tPrime_helper_params']['tPrime_path'], 'runit.bat')
    elif sys.platform.startswith('linux'):
        exe_path = os.path.join(args['tPrime_helper_params']['tPrime_path'], 'runit.sh')
    else:
        print('unknown system, cannot run TPrime')   
        
    # Print out command for help with debugging
    tcmd = exe_path + ' -syncperiod=' + repr(sync_period) + \
        ' -tostream=' + toStream_path
        
    if offset_exists:
        tcmd = tcmd + ' -offsets=' + offset_path

    for i, fp in enumerate(from_list):
        tcmd = tcmd + ' -fromstream=' + repr(i) + ',' + fp

    for i, ep in enumerate(events_list):
        tcmd = tcmd + ' -events=' + repr(from_stream_index[i]) + ',' + ep + ',' + out_list[i]

    # write out file to record the TPrime command for a record
    bat_path = os.path.join(run_directory, run_name + '_TPrime_cmd.txt')
    with open(bat_path, 'w') as batfile:
        batfile.write(tcmd)

        
    # make the TPrime call
    subprocess.Popen(tcmd,shell='False').wait()

    # convert output files were text, convert to npy
    if not bNPY:
        for op in out_list:
            spike_times_sec_to_npy(op)
            
            
    # if the psth extract string is not null
    
    extract_str = args['tPrime_helper_params']['psth_ex_str']
    if len(extract_str) > 0:       
        create_PSTH_events( all_list, out_list, sort_out_list, extract_str, \
                           args['tPrime_helper_params']['sort_out_tag'] )

    execution_time = time.time() - start

    print('total time: ' + str(np.around(execution_time, 2)) + ' seconds')

    return {"execution_time": execution_time}  # output manifest


def spike_times_npy_to_sec(sp_fullPath, sample_rate, bNPY):
    # convert spike_times.npy to text of times in sec
    # return path to the new file. Can take sample_rate as a
    # parameter, or set to 0 to read from param file

    # get file name and create path to new file
    sp_path, sp_fileName = os.path.split(sp_fullPath)
    baseName, bExt = os.path.splitext(sp_fileName)
    if bNPY:
        new_fileName = baseName + '_sec.npy'
    else:
        new_fileName = baseName + '_sec.txt'
        
    new_fullPath = os.path.join(sp_path, new_fileName)

    # load spike_times.npy; returns numpy array (Nspike,) as uint64
    spike_times = np.load(sp_fullPath)

    if sample_rate == 0:
        # get sample rate from params.py file, assuming sp_path is a full set
        # of phy output
        with open(os.path.join(sp_path, 'params.py'), 'r') as f:
            currLine = f.readline()
            while currLine != '':  # The EOF char is an empty string
                if 'sample_rate' in currLine:
                    sample_rate = float(currLine.split('=')[1])
                    print(f'sample_rate read from params.py: {sample_rate:.10f}')
                currLine = f.readline()

            if sample_rate == 0:
                print('failed to read in sample rate\n')
                sample_rate = 30000

    spike_times_sec = spike_times/sample_rate   # spike_times_sec dtype = float

    if bNPY:
        # write out npy file
        np.save(new_fullPath, spike_times_sec)
    else:
        # write out single column text file
        nSpike = len(spike_times_sec)
        with open(new_fullPath, 'w') as outfile:
            for i in range(0, nSpike-1):
                outfile.write(f'{spike_times_sec[i]:.6f}\n')
            outfile.write(f'{spike_times_sec[nSpike-1]:.6f}')

    return new_fullPath


def spike_times_sec_to_npy(spa_fullPath):
    # convert a text file of spike times in seconds to an npy file of
    # python floats, with shape (Nspike,)
    # spa => spikes, adjusted

    # get file name and create path to new file
    spa_path, spa_fileName = os.path.split(spa_fullPath)
    baseName, bExt = os.path.splitext(spa_fileName)
    new_fileName = baseName + '.npy'
    new_fullPath = os.path.join(spa_path, new_fileName)

    # count lines in file to size array
    lineCount = 0
    with open(os.path.join(spa_fullPath), 'r') as f:
        for line in f:
            lineCount = lineCount + 1

    times = np.zeros((lineCount,), dtype='float')
    # read in text file, single column of floats
    i = 0
    with open(os.path.join(spa_fullPath), 'r') as f:
        for line in f:
            times[i] = float(line)
            i = i + 1

    np.save(new_fullPath, times)
    
def parse_stream(stream_str):
    # stream_str is the stream name followed by the index (e.g. imec3)
    # stream type indicies (js): 0 = ni, 1 = obx, 2 = imec    
    if stream_str.find('ni') > -1:
        js = 0
        ip = 0 # never more than one ni stream
    elif stream_str.find('obx') > -1:
        js = 1
        temp = stream_str.partition('obx')[2]
        ip = int(temp.partition('_')[0])
    elif stream_str.find('imec') > -1:
        js = 2
        ip = int(stream_str.partition('imec')[2])
    else:
        js = -1
        ip = -1
    return js, ip  
        
    
    
def parse_catgt_fyi(fyi_path, toStream_id):

    # read fyi file, build array of paths to sync files 
    # toStream_id is the tuple (js,ip) for the toStream
    # return 
    #   toStream_path
    #   from_list = list of paths to sync edge files for each from stream
    #   event_list = list of event edges in the fyi file
    #   fron_stream_index = list of from stream indicies for fyi events files
    #   out_list = paths for output for fyi events files
    #          
    
    # create empty lists to fill in
    from_list = []
    from_list_ids = []
    events_list = []
    from_stream_index = []
    out_list = []
    all_list = []
   
    with open(fyi_path, 'r') as reader:
        line = reader.readline()
        while line != '':  # The EOF char is an empty string
            
            stream_str, eq, curr_path = line.partition("=")
            stream_str = stream_str.partition('_')[2]
            
            if line.find('sync') == 0:                
                # this is the path to a file of sync edges
                curr_path = curr_path[0:len(curr_path)-1] # trim off cr on end
                all_list.append(curr_path)
                js, ip = parse_stream(stream_str)
                if (js,ip) == toStream_id:
                    toStream_path = curr_path
                else:
                    from_list_ids.append((js,ip))
                    from_list.append(curr_path)
  
            elif line.find('times') == 0:
                # this is the path to a file of event edges
                # if it is not from the toStream, add to events, and out_list
                curr_path = curr_path[0:len(curr_path)-1] # trim off cr on end
                all_list.append(curr_path)
                js, ip = parse_stream(stream_str)
                if (js,ip) != toStream_id:
                    events_list.append(curr_path)
                    from_stream_index.append(from_list_ids.index((js,ip)))
                    cp = Path(curr_path)
                    curr_output_name = cp.stem + '.adj' + cp.suffix            
                    out_list.append(os.path.join(cp.parent, curr_output_name))
            line = reader.readline()   
         
    return toStream_path, from_list, from_list_ids, events_list, from_stream_index, out_list, all_list


def create_prbDir_list(run_directory, prb_dir_prefix):
    # Make a list of all probe directories in run_directory
    # used by create_PSTH_events to copy teh events file to each probe directory
    
    prb_fld_wild = prb_dir_prefix + '*'
    prbDir_list = list()
    for pDir in os.listdir(run_directory):
        if fnmatch.fnmatch(pDir, prb_fld_wild):
            prbDir_list.append(os.path.join(run_directory, pDir))

    return prbDir_list


def create_PSTH_events( all_list, out_list, prbDir_list, extract_str, sort_name ):
    # For a user specified event set, create a csv file and copy to all kilosort
    # output folders. Get a search string based on the extract_str provided by the user
    ex_type, stream_index, prb_index, ex_name_str = catGT_ex_params_from_str(extract_str)
    search_pat_adj = '*' + ex_name_str + '.adj.txt'
    search_pat_orig = '*' + ex_name_str + '.txt'
      
    # check out_list for adjusted file
    found = False
    n_out = len(out_list)
    n = 0
    while not found and n < n_out:       
        eventPath = out_list[n]
        print(repr(n) + ': ' + eventPath)
        if fnmatch.fnmatch(eventPath, search_pat_adj):
            found = True
        else:
            n = n + 1
    # if there isn't an adjusted file, get the path from list of all
    n = 0
    n_all = len(all_list)
    while not found and n < n_all:
        eventPath = all_list[n]
        print(repr(n) + ': ' + eventPath)
        if fnmatch.fnmatch(eventPath, search_pat_orig):
            found = True
        else:
            n = n + 1 
            
    print(eventPath)
    
    if found:
    # the CatGT extracted edge files are a single column with </n>
    # event viewer needs .csv
        edgeTimes = np.zeros((0), dtype='float')
        with open(eventPath, 'r') as inFile:
            line = inFile.readline()
            while line != '':  # The EOF char is an empty string
                currEdge = float(line)
                edgeTimes = np.append(edgeTimes, currEdge)
                line = inFile.readline()

    # The output should be saved with the phy output, where the event viewer
    # plugin can read it.   
        for pDir in prbDir_list:
            event_path = os.path.join(pDir, 'events.csv')
            nEvent = len(edgeTimes)
            with open(event_path, 'w') as outfile:
                if nEvent > 0:
                    for i in range(0, nEvent-1):
                        outfile.write(f'{edgeTimes[i]:.6f},')
                    outfile.write(f'{edgeTimes[nEvent-1]:.6f}')

    return
        

def main():

    from ._schemas import InputParameters, OutputParameters

    """Main entry point:"""
    mod = ArgSchemaParser(schema_type=InputParameters,
                          output_schema_type=OutputParameters)

    output = call_TPrime(mod.args)

    output.update({"input_parameters": mod.args})
    if "output_json" in mod.args:
        mod.output(output, indent=2)
    else:
        print(mod.get_output_json(output))


if __name__ == "__main__":
    main()
