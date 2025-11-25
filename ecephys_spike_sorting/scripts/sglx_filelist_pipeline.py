import os
import sys
from helpers import SpikeGLX_utils
from helpers import run_one_probe
from create_input_json import createInputJson


# script to run sorting and postprocessing moduels on a list of recordings
# the SpikeGLX metadata must be in the directory with the binary.
# phy output is directed to a newly made directory next to the input binary.
# the primary use for this script is re-running sorting on data that has been
# preprocessed by CatGT


# -------------------------------
# -------------------------------
# User input -- Edit this section
# -------------------------------
# -------------------------------
# Full path to log file, including file name
# If this file exists, new run data is appended to it
logFullPath = r'D:\SC048_out\pipeline_log.csv'

# ks_ver  sets up the output tag and threshold values.
# To run a specific MATLAB KS, make sure to set up the kilosort_repository in 
# create_input_json; in the module list, call 'kilosort helper'
# To run KS4, use an anaconda install and call 'ks4_helper'

ks_ver = '4'  # needs to be one of: '2.0', '2.5', '3.0', or '4'
ksTag_dict = {'2.0':'ks2', '2.5':'ks25', '3.0':'ks3', '4':'ks4'}
ks_output_tag = ksTag_dict[ks_ver]

# brain region specific params
# can add a new brain region by adding the key and value for each param
# can add new parameters -- any that are taken by create_input_json --
# by adding a new dictionary with entries for each region and setting the 
# according to the new dictionary in the loop to that created json files.
#

refPerMS_dict = {'default': 2.0, 'cortex': 2.0, 'medulla': 1.5, 'thalamus': 1.0, 'striatum': 2.0, 'midbrain':1.5}

# threhold values appropriate for KS2, KS2.5
ksTh2_dict = {'default':'[10,4]', 'cortex':'[10,4]', 'medulla':'[10,4]', 'thalamus':'[10,4]'}
# threshold values appropriate for KS3.0
ksTh3_dict = {'default':'[9,9]', 'cortex':'[9,9]', 'medulla':'[9,9]', 'thalamus':'[9,9]'}
# threshold values appropriate for KS4.0
ksTh4_dict = {'default':'[8,9]', 'cortex':'[8,9]', 'medulla':'[8,9]', 'thalamus':'[8,9]', 'striatum':'[8,9]', 'midbrain':'[8,9]'}

if ks_ver == '2.0' or ks_ver == '2.5':
    ksTh_dict = ksTh2_dict
elif ks_ver == '3.0':    
    ksTh_dict = ksTh3_dict
elif ks_ver == '4':
    ksTh_dict = ksTh4_dict
else:
    print('unknown version of ks, exiting.')
    sys.exit()

# -----------
# Input data
# -----------


# for each recording, specfiy a full path the the binary and a brain region

recording_specs = [		
    
                [r'D:\SC048_out\catgt_SC048_122920_ex_g0\SC048_122920_ex_g0_imec2\SC048_122920_ex_g0_tcat.imec2.ap.bin',['cortex']]

]



# --------------------------
# KS2, KS2.5, KS3 parameters
# --------------------------
# parameters that will be constant for all recordings
# Template ekmplate radius and whitening, which are specified in um, will be 
# translated into sites using the probe geometry.
ks_remDup = 0       # used by KS2, 2.5, 3
ks_saveRez = 1      # used by KS2, 2.5, 3
ks_copy_fproc = 0   # used by 2.5, 3, to save drift corrected binary
ks_templateRadius_um = 163    # used by KS2, 2.5, 3
ks_whiteningRadius_um = 163   # used by KS2, 2,5 2.5, 3
ks_minfr_goodchannels = 0.1   # used by KS2, 2.5, 3; set to 0 for KS2.5 and 3

# -------------------------------
# KS2, KS2.5, KS3, KS4 parameters
# -------------------------------
ks_CAR = 0          # CAR already done in catGT
ks_nblocks = 6      # for KS2.5 KS3, and KS4; 1 for rigid registration in drift correction, 
                    # higher numbers to allow different drift for different 'blocks' of the probe

# -------------------------------------------------------
# KS4 specific parameters -- these are the default values
# -------------------------------------------------------
ks4_duplicate_spike_ms = 0.25
ks4_min_template_size_um = 10

# If running KS20_for_preprocessed_data:
# (https://github.com/jenniferColonell/KS20_for_preprocessed_data)
# can skip filtering with the doFilter parameter.
# Useful for speed when data has been filtered with CatGT.
# This parameter is not implemented in standard versions of kilosort.
ks_doFilter = 0

# ----------------------
# C_Waves snr radius, um
# ----------------------
c_Waves_snr_um = 160



# ---------------
# Modules List
# ---------------
# List of modules to run per probe; CatGT and TPrime are called once for each run.
modules = [
			'ks4_helper',
            'kilosort_postprocessing',
            #'noise_templates',
            'mean_waveforms',
            'quality_metrics'
			]

json_directory = r'D:\SC048_out\json_files'


# -----------------------
# -----------------------
# End of user input
# -----------------------
# -----------------------
if ks_ver in ['2.0','2.5','3.0'] and 'kilosort_helper' not in modules and 'ks4_helper' in modules:
    print('For MATLAB versions of KS, run kilosort_helper module')
    sys.exit()
if ks_ver == '4' and 'ks4_helper' not in modules and 'kilosort_helper' in modules:
    print('For kilsort 4, run ks4_helper module')
    sys.exit()

# delete existing C_waves.log
try:
    os.remove('C_Waves.log')
except OSError:
    pass


# first loop over recording specs to create output directories and
#  make json file for each

# initialize lists    
module_input_json = []
module_output_json = []
session_id = []
data_directory = []
output_dir = []

kilosort_output_parent = r'D:\SC048_out\filelist_output'
 
for i, spec in enumerate(recording_specs):
    
    path = spec[0]
    npx_directory = os.path.dirname(path)
    fname, fextension = os.path.splitext(path)
    input_meta_fullpath = os.path.join(npx_directory, (fname + '.meta'))
    print(input_meta_fullpath)
    binName = os.path.basename(path)
    baseName = SpikeGLX_utils.ParseTcatName(binName)
    prbStr = SpikeGLX_utils.GetProbeStr(binName)   # returns empty string for 3A


    # Asuume that ouput parent = input
    
    kilosort_output_run = os.path.join(kilosort_output_parent,f'{baseName}_imec{prbStr}')
    # print(kilosort_output_run)
    if not os.path.exists(kilosort_output_run):
       os.mkdir(kilosort_output_run)
        
        
    # output subdirectory
    outputName = 'imec' + prbStr + '_' + ks_output_tag
    
    kilosort_output_dir = os.path.join(kilosort_output_run, outputName)
    
    session_id.append(f'{baseName}_imec{prbStr}') 
    
    module_input_json.append(os.path.join(json_directory, session_id[i] + '-input.json'))
    
    data_directory.append(npx_directory)
    continuous_file = os.path.join(data_directory[i], binName)
 
 

    # kilosort_postprocessing and noise_templates moduules alter the files
    # that are input to phy. If using these modules, keep a copy of the
    # original phy output
    if ('kilosort_postprocessing' in modules) or('noise_templates' in modules):
        ks_make_copy = True
    else:
        ks_make_copy = False

    
    # get region specific parameters
    ks_Th = ksTh_dict.get(spec[1][0])
    refPerMS = refPerMS_dict.get(spec[1][0])
   

    info = createInputJson(module_input_json[i], npx_directory=npx_directory,
                                   continuous_file = continuous_file,
                                   spikeGLX_data = True,
                                   input_meta_path = input_meta_fullpath,
                                   kilosort_output_directory = kilosort_output_dir,
                                   ks_make_copy = ks_make_copy,
                                   noise_template_use_rf = False,
                                   catGT_run_name = session_id[i], 
                                   ks_ver = ks_ver,
                                   ks_remDup = ks_remDup,                   
                                   ks_finalSplits = 1,
                                   ks_labelGood = 1,
                                   ks_saveRez = ks_saveRez,
                                   ks_copy_fproc = ks_copy_fproc,
                                   ks_minfr_goodchannels = ks_minfr_goodchannels,                  
                                   ks_whiteningRadius_um = ks_whiteningRadius_um,
                                   ks_doFilter = ks_doFilter,
                                   ks_Th = ks_Th,
                                   ks_CSBseed = 1,
                                   ks_LTseed = 1,
                                   ks_nblocks = ks_nblocks,
                                   ks_templateRadius_um = ks_templateRadius_um,
                                   extracted_data_directory = npx_directory,
                                   c_Waves_snr_um = c_Waves_snr_um,                               
                                   qm_isi_thresh = refPerMS/1000,
                                   ks4_duplicate_spike_ms = ks4_duplicate_spike_ms,
                                   ks4_min_template_size_um = ks4_min_template_size_um,
                                   include_pc_metrics = True
                                   
                                   )   

    # copy json file to data directory as record of the input parameters 
       
        
# loop over files again for processing. 
# not running CatGT; set params accorindly
run_CatGT = False
catGT_input_json = ''
catGT_output_json = ''

for i, spec in enumerate(recording_specs):
       
    run_one_probe.runOne( session_id[i],
             json_directory,
             data_directory[i],
             run_CatGT,
             catGT_input_json,
             catGT_output_json,
             modules,
             module_input_json[i],
             logFullPath )
                 
        
