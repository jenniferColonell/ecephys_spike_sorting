# -*- coding: utf-8 -*-
import numpy as np
import fnmatch
import os
import sys
from pathlib import Path


import ecephys_spike_sorting.common.SGLXMetaToCoords as SGLXMeta

def GetFirstTrialPath(catGT_run_name, gate_string, trigger_string, probe_string ):
    prb_list = ParseProbeStr(probe_string)
    run_folder = catGT_run_name + '_g' + gate_string
    prb_folder =  run_folder + '_imec' + prb_list[0]
    first_trig, last_trig = ParseTrigStr(trigger_string, prb_folder)
    filename = catGT_run_name + '_g' + gate_string + '_t' + \
                 str(first_trig) + '.imec' + prb_list[0] + '.ap.bin'
    firstTrialPath = os.path.join(run_folder, prb_folder, filename) 

    return firstTrialPath

def GetTrialRange(prb, gate, prb_folder):
    tFiles = os.listdir(prb_folder)
    minIndex =  sys.maxsize
    maxIndex = 0
    searchStr = '_g' + gate + '_t'
    print(searchStr)
    for tName in tFiles:
        if (fnmatch.fnmatch(tName,'*.bin')):
            gPos = tName.find(searchStr)
            tStart = gPos + len(searchStr)
            tEnd = tName.find('.', tStart)
            
            if gPos > 0 and tEnd > 0:
                try:
                    tInd = int(tName[tStart:tEnd])                    
                except ValueError:
                    print(tName[tStart:tEnd])
                    print('Error parsing trials for probe folder: ' + prb_folder + '\n')
                    return -1, -1
            else:
                print('Error parsing trials for probe folder: ' + prb_folder + '\n')
                return -1, -1
            
            if tInd > maxIndex:
                maxIndex = tInd
            if tInd < minIndex:
                minIndex = tInd
                
    return minIndex, maxIndex

    
def EphysParams(metaFullPath):
    # get ephys params from metadata at meta full path    
    # read metadata
    
    #first create Path object from string
    metaPath = Path(metaFullPath)
    meta = SGLXMeta.readMeta(metaPath)
    
    if 'imDatPrb_type' in meta:
        pType = (meta['imDatPrb_type'])
        if pType =='0':
            probe_type = 'NP1'
        else:
            probe_type = 'NP' + pType
    else:
        probe_type = '3A'    #3A probe
    
    sample_rate = float(meta['imSampRate'])    
    
    num_channels = int(meta['nSavedChans'])
    
    uVPerBit = Chan0_uVPerBit(meta, probe_type)
    
    # get number of AP and SY channels in file
    # used in waveform metrics grab only AP channels for calculations
    # Note that with quad base SY can be larger than 1
    nAP, nLF, nSY = SGLXMeta.ChannelCountsIM(meta)
    
    if 'snsGeomMap' in meta:
        useGeom = True
    else:
        useGeom = False
    
        
    # read shank map to get disabled (reference) channels
    ref_channels = GetDisabledChan(meta, useGeom)
    
    xCoord, yCoord, shankInd, connected = SGLXMeta.MetaToCoords(metaPath,-1)
    sh, sh_counts = np.unique(shankInd, return_counts=True)
    # get vpitch, hpitch, nColumn from shank with the largest numbers of sites
    sh_mode = sh[np.argmax(sh_counts)]
    sh_mode_sites = np.squeeze( np.argwhere(shankInd == sh_mode))
    x_sh = xCoord[sh_mode_sites]
    y_sh = yCoord[sh_mode_sites]
    nColumn = len(np.unique(x_sh))    
    
    # NP 1.0 staggered patterns have been historically as 2 columns. These can
    # be tricky to identify from the coorindates alone AND there are some type 0 
    # probes that are NOT staggered.
    vpitch = 0
    hpitch = 0
    stag_types = ['3A', 'NP1', 'NP1200', 'NP1210', 'NP1020', 'NP1021', 'NP1030', 'NP1031']
    if (probe_type in stag_types) and (nColumn == 4):
        nColumn = 2
        vpitch = 20  # happens to be true for all staggered configurations
        if probe_type in ['NP1020', 'NP1021', 'NP1030', 'NP1031']:
            hpitch = 87
        else:
            hpitch = 32
    else:
        # these are linear or square pattern probes        
        xval, xval_counts = np.unique(x_sh, return_counts=True)        
        xval_mode = xval[np.argmax(xval_counts)]
        xval_mode_sites = np.squeeze( np.argwhere(x_sh == xval_mode))        
        yval_col = y_sh[xval_mode_sites]        
        if len(yval_col) > 1:
            ydiff = np.diff(np.sort(yval_col))
            ydiffval, ydiff_counts = np.unique(ydiff, return_counts=True)
            vpitch = ydiff[np.argmax(ydiff_counts)]
        if len(xval) > 1:
            xdiff = np.diff(np.sort(xval))
            xdiffval, xdiff_counts = np.unique(xdiff, return_counts=True)
            hpitch = xdiff[np.argmax(xdiff_counts)]
    
    return(probe_type, sample_rate, num_channels, ref_channels, uVPerBit, vpitch, hpitch, nColumn, nAP, nSY, useGeom)


# Return gain for imec channels.
# Index into these with the original (acquired) channel IDs.
#
def Chan0_uVPerBit(meta, probe_type):
    # Returns uVPerBit conversion factor for channel 0
    # If all channels have the same gain (usually set that way for 
    # 3A and NP1 probes; always true for NP2 probes), can use
    # this value for all channels.
    
    # first check if metadata includes the imChan0apGain key
    if 'imChan0apGain' in meta:
        APgain = float(meta['imChan0apGain'])
        voltage_range = float(meta['imAiRangeMax']) - float(meta['imAiRangeMin'])
        maxInt = float(meta['imMaxInt'])
        uVPerBit = (1e6)*(voltage_range/APgain)/(2*maxInt)
        
    else:     
        imroList = meta['imroTbl'].split(sep=')')
        # One entry for each channel plus header entry,
        # plus a final empty entry following the last ')'
        # channel zero is the 2nd element in the list
    
        if probe_type == 'NP21' or probe_type == 'NP24':
            # NP 2.0; APGain = 80 for all channels
            # voltage range = 1V
            # 14 bit ADC
            uVPerBit = (1e6)*(1.0/80)/pow(2,14)
        elif probe_type == 'NP1110':
            # UHD2 with switches, special imro table with gain in header            
            currList = imroList[0].split(sep=',')
            APgain = float(currList[3])
            uVPerBit = (1e6)*(1.2/APgain)/pow(2,10)
        else:
            # 3A, 3B1, 3B2 (NP 1.0), or other NP 1.0-like probes
            # voltage range = 1.2V
            # 10 bit ADC
            currList = imroList[1].split(sep=' ')   # 2nd element in list, skipping header
            APgain = float(currList[3])
            uVPerBit = (1e6)*(1.2/APgain)/pow(2,10)
        
    return(uVPerBit)


def GetDisabledChan(meta, useGeom):
    
    chanCountList = meta['snsApLfSy'].split(sep=',')
    AP = int(chanCountList[0])
    
    if useGeom is True:
        useMap = meta['snsGeomMap'].split(sep=')')
    else:
        useMap = meta['snsShankMap'].split(sep=')')
        
    # loop over known number of AP channels to avoid problems with
    # extra entries in older data
    connected = np.zeros((AP,)) 
    for i in range(AP):
        # get parameter list from this entry, skipping first header entry
        currEntry = useMap[i+1]
        currList = currEntry.split(sep=':')
        connected[i] = int(currList[3]) 
        
    disabled_chan = np.where(connected==0)[0].tolist()
    
    return disabled_chan


def ParseProbeStr(probe_string):
    
    str_list = probe_string.split(',')
    prb_list = []
    for substr in str_list:
        if (substr.find(':') > 0):
            # split at colon
            subsplit = substr.split(':')
            for i in range( int(subsplit[0]), int(subsplit[1]) + 1):
                prb_list.append(str(i))
        else:
            # just append this string
            prb_list.append(substr)

    return prb_list

def ParseTrigStr(trigger_string, prb, gate, prb_folder):
    
    str_list = trigger_string.split(',')
    first_trig_str = str_list[0]
    last_trig_str = str_list[1]
    
    if last_trig_str.find('end') >= 0 or first_trig_str.find('start') >= 0 :
        # get the full range from the directory
        minInd, maxInd = GetTrialRange(prb, gate, prb_folder)

    if first_trig_str.find('start') >= 0:
        first_trig = minInd
    else:
        first_trig = int(first_trig_str)
    
    if last_trig_str.find('end') >= 0:
        last_trig = maxInd
    else:
        last_trig = int(last_trig_str)
        
    # trig_array =  np.arange(first_trig, last_trig+1)

    return first_trig, last_trig


def ParseGateStr(gate_string):
    str_list = gate_string.split(',')
    first_gate = int(str_list[0])
    if len(str_list) > 1:
        last_gate = int(str_list[1])
    else:
        last_gate = first_gate
    return first_gate, last_gate

def ParseTcatName(tcat_name):    
    tcat_pos = tcat_name.find('tcat',0)
    baseName = tcat_name[0:tcat_pos-1]  #subtrace 1 from tcat pos to remove _
    return baseName

def GetProbeStr(tcat_name):
    tcat_pos = tcat_name.find('tcat',0)
    ap_pos = tcat_name.find('.ap',0)   
    imStr = tcat_name[tcat_pos+5:ap_pos]
    if len(imStr) == 4:
        prbStr = ''      # 3A data, no probe index
    else:
        prbStr = imStr[4:len(imStr)]
    return prbStr


def ParseCatGTLog(logPath, run_name, gate_string, prb_list):

    gfix_str = run_name + '_' + gate_string + ' Gfix'

    num_probe = len(prb_list)
    gfix_edits = np.zeros(num_probe, dtype='float64')

    gfound = np.zeros(num_probe)
    pfound = list()             # list of strings of probes found
    nfound = 0
    log_fullpath = logPath.replace('\\', '/') + "/CatGT.log"

    with open(log_fullpath, 'r') as reader:
        line = reader.readline()
        while line != '' and nfound < num_probe:  # The EOF char is an empty string
            gstart = line.find( gfix_str )
            if gstart  > -1:      
                # parse this line to get probe string and corrections/sec
                line_end = len(line)
                gsub = line[gstart:line_end]
                datArr = gsub.split()       
                pfound.append(datArr[3])
                gfound[nfound] = float(datArr[5])
                nfound = nfound + 1
            line = reader.readline()   
    
    # order the returned gfix_edits matching the probe order specified 
    # in prb_list
    for i in range(0,len(prb_list)):
        if prb_list[i] in pfound:
            gfix_edits[i] = gfound[pfound.index(prb_list[i])]
     
    return gfix_edits

def CreateAuxTimeEvents(catGT_run_name, first_gate_string, catGT_dest, catGT_out_tag, stream="ni"):
    # Creates a set events are simply the times for the collected data 
    # in any auxiliary (non-neural) stream, specified by stream.
    # 
    # CatGT 1.9 and later always creates an output NI metadata file
    # 
 
    output_folder = f'{catGT_out_tag}_{catGT_run_name}_g{first_gate_string}'
         
    if stream.find('ni') > -1:
        auxMeta_filename = f'{catGT_run_name}_g{first_gate_string}_tcat.nidq.meta'
        auxMeta_path = Path(os.path.join(catGT_dest, output_folder, auxMeta_filename))       
        meta = SGLXMeta.readMeta(auxMeta_path)
        sample_rate = float(meta['niSampRate'])
    elif stream.find('ob') > -1:
        auxMeta_filename = f'{catGT_run_name}_g{first_gate_string}_tcat.{stream}.obx.meta'
        auxMeta_path = Path(os.path.join(catGT_dest, output_folder, auxMeta_filename))       
        meta = SGLXMeta.readMeta(auxMeta_path)
        sample_rate = float(meta['obSampRate'])
    else:
        print('error in CreateAuxTimeEvents, unknow stream type')
        return
    
    num_channels = int(meta['nSavedChans'])
    nSamp = int(meta['fileSizeBytes'])/(num_channels * 2)
    aux_times = np.arange(nSamp)/sample_rate
    
    # save aux_times in output folder to be an event file

    out_name = f'{catGT_run_name}_g{first_gate_string}_tcat.{stream}.times.npy'
    out_path = os.path.join(catGT_dest, output_folder, out_name)
    np.save(out_path,aux_times)
    
    # check for presence of an all_fyi file, indicating run with catgt 3.0 or later in pipeline
    # supercat output will be named _fyi.txt 
    fyi_all_path = Path(os.path.join(catGT_dest, output_folder, catGT_run_name + '_g' + first_gate_string + '_all_fyi.txt'))
    fyi_single_path = Path(os.path.join(catGT_dest, output_folder, catGT_run_name + '_g' + first_gate_string + '_fyi.txt'))
    fyi_exists = False
    if Path(fyi_all_path).is_file():
        fyi_exists = True
        fyi_path = fyi_all_path
    elif Path(fyi_single_path).is_file():
        fyi_exists = True
        fyi_path = fyi_single_path

    if fyi_exists:
        # append a line for the newly created times file        
        file_fyi = open(fyi_path, "a")  # append mode
        file_fyi.write(f'times_{stream}_N=' + out_path + '\n')
        file_fyi.close()

    return

def Chans2PrintStr(chan_bool):
    # return printer string (runs of sequential channels given by a:b) for 
    # an array of channel indicies
    # code adapted from Bill Karsh
    # create a logical array with maxChan elements
    chan_ind = np.where(chan_bool)
    nb = len(chan_bool)
    min_ind = np.min(chan_ind)

    out_str = repr(min_ind)
    rn = False
    
    # starting from the first channel
    for ib in range(min_ind+1, nb):
        if not chan_bool[ib]:
             # if we were in a run, the last index was the end of it
            if rn:
                out_str = out_str + ':' + repr(ib-1)
                rn = False
        else:
            if not chan_bool[ib-1]:
                # if the last index was not included start a new entry
                out_str = out_str + ',' + repr(ib)
            else:
                # if the last index was included, we have a run
                rn = True
    
    # after looping over all indices, add last entry if in a run
    if rn:
        out_str = out_str + ':' + repr(nb-1)            
    
    return out_str
    
def CreateShankSaveString(metaFullPath):
    # superseded by CreateSepShanksString with CatGT version 4.4
    # first create Path object from string
    metaPath = Path(metaFullPath)
    metaStem = metaPath.stem
    imec_pos = metaStem.rfind('imec',0)
    ap_pos = metaStem.rfind('.ap',0) 
    if imec_pos > 0 and ap_pos > 0:
        # good metadata name
        imStr = metaStem[imec_pos:ap_pos]
        if len(imStr) == 4:
            prb_ind = 0      # 3A data, no probe index
        else:
            prb_ind = int(imStr[4:len(imStr)])
    else:
       print('Incorrect metadata filenmae.') 
       nShank = 0
       save_str = ''
       out_ind_list = []
       return nShank, save_str, out_ind_list
    
    # check for a sync string in the input
    meta = SGLXMeta.readMeta(metaPath)
    AP, LF, SY = SGLXMeta.ChannelCountsIM(meta)
    orig_chan_ind = SGLXMeta.OriginalChans(meta)    
    ap_chan_ind = orig_chan_ind[0:AP]
    sy_ind = orig_chan_ind[len(orig_chan_ind)-SY:len(orig_chan_ind)]
    print(repr(sy_ind))

    xCoord, yCoord, shankInd, connected = SGLXMeta.MetaToCoords(metaPath,-1)    
    sh, sh_counts = np.unique(shankInd, return_counts=True)
    
    save_str = ''
    out_ind_list = []
    nShank = len(sh)
    for i in range(nShank):
        # get channels on this shank
        sh_chan_bool = (shankInd == sh[i])
        all_chan_bool = np.ndarray((np.max(orig_chan_ind)+1,),dtype=bool)
        all_chan_bool[:] = False
        all_chan_bool[sy_ind] = True
        all_chan_bool[ap_chan_ind[sh_chan_bool]] = True
        sh_chan_str = Chans2PrintStr(all_chan_bool)        
        out_ind_str = str(int(1000 + 10*prb_ind + sh[i]))
        out_ind_list.append(out_ind_str)
        save_str = save_str + ' -save=2,' + repr(prb_ind) + ',' + out_ind_str + ',' + sh_chan_str
        
    return nShank, save_str, out_ind_list

def CreateSepShanksString(metaFullPath):
    # create sepShanks string for CatGT 44 and later
    # The sepShanks option automatically generates a save 
    # string with the channels for each shank.
    # For the 2020 probe (Quad base) the SYNC signal for each shank
    # is saved with the neural data. For other multishank probes, the 
    # common SYNC signal is saved with each shank.
    # first create Path object from string
    metaPath = Path(metaFullPath)
    metaStem = metaPath.stem
    imec_pos = metaStem.rfind('imec',0)
    ap_pos = metaStem.rfind('.ap',0) 
    if imec_pos > 0 and ap_pos > 0:
        # good metadata name
        imStr = metaStem[imec_pos:ap_pos]
        if len(imStr) == 4:
            prb_ind = 0      # 3A data, no probe index
        else:
            prb_ind = int(imStr[4:len(imStr)])
    else:
       print('Incorrect metadata filenmae.') 
       sh = []
       sepShanks_str = ''
       out_ind_list = []
       
    nShank, shankInd = shank_info(metaPath)
    sh = np.unique(shankInd)
    sepShanks_str = f'-sepShanks={prb_ind}'
    out_ind_list = []
    
    for i in range(nShank):  
        if np.isin(i,sh):
            out_ind_str = str(int(1000 + 10*prb_ind + i))
            out_ind_list.append(out_ind_str)
        else:
            out_ind_str = '-1'
        sepShanks_str = sepShanks_str + ','  + out_ind_str

        
           
    return len(sh), sepShanks_str, out_ind_list


def shank_info(metaFullPath):
    # re-reads the metadata, which is inefficient.
    # However, keeps the calling script simple
        
    metaPath = Path(metaFullPath)
    meta = SGLXMeta.readMeta(metaPath)
    
    if 'snsGeomMap' in meta:
        [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = SGLXMeta.geomMapToGeom(meta)
    else:   
        [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = SGLXMeta.shankMapToGeom(meta)
    
    return nShank, shankInd

def rec_length(metaFullPath):
    # re-reads the metadata, which is inefficient.
    # However, keeps the calling script simple (returns only the time)
    
    # first create Path object from string
    metaPath = Path(metaFullPath)
    meta = SGLXMeta.readMeta(metaPath)
    
    # use filesize, saved channels, and rate,
    # just in case this is a user-edited file that didn't get  fileTimeSecs set correctly
    nSamp = int(meta['fileSizeBytes'])/(2*int(meta['nSavedChans']))
    recTimeSec = nSamp/float(meta['imSampRate'])
    return recTimeSec