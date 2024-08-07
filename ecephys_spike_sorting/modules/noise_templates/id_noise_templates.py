import numpy as np

from scipy.signal import correlate, find_peaks, cwt, ricker
from sklearn.ensemble import RandomForestClassifier

from scipy.interpolate import griddata
from scipy.ndimage.filters import gaussian_filter1d

from ...common.utils import printProgressBar

import multiprocessing
from functools import partial

import pickle

def id_noise_templates_rf(spike_times, spike_clusters, cluster_ids, templates, params):

    """
    Uses a random forest classifier to identify noise units based on waveform shape

    Inputs:
    -------
    spike_times : spike times (in seconds)
    spike_clusters : cluster IDs for each spike time []
    cluster_ids : all unique cluster ids
    templates : template for each unit output by Kilosort

    Outputs:
    -------
    cluster_ids : same as input
    is_noise : boolean array, True at index of noise templates

    Parameters:
    ----------
    'classifier_path' : path to pickled classifier object

    """
    
    # #############################################

    classifier_path = params['classifier_path']

    # #############################################

    classifier = pickle.load(open(classifier_path, 'rb'))

    feature_matrix = np.zeros((cluster_ids.size, 61, 32))

    peak_channels = np.squeeze(np.argmax(np.max(templates,1) - np.min(templates,1),1))

    for idx, unit in enumerate(cluster_ids):
        
        peak_channel = peak_channels[unit]

        min_chan = np.max([0,peak_channel-16])
        if min_chan == 0:
            max_chan = 32
        else:
            max_chan = np.min([templates.shape[2], peak_channel+16])
            if max_chan == templates.shape[2]:
                min_chan = max_chan - 32

        sub_template = templates[unit, :, min_chan:max_chan]

        feature_matrix[idx,:,:] = sub_template

    feature_matrix = np.reshape(feature_matrix[:,:,:], (feature_matrix.shape[0], feature_matrix.shape[1] * feature_matrix.shape[2]), 2)
    feature_matrix = feature_matrix[:,::4]

    is_noise = classifier.predict(feature_matrix)
    is_noise = is_noise.astype('bool')

    return cluster_ids, is_noise
    


def id_noise_templates(cluster_ids, templates, channel_pos, params):

    """
    Uses a set of heuristics to identify noise units based on waveform shape

    Inputs:
    -------
    cluster_ids : all unique cluster ids
    templates : template for each unit output by Kilosort
    channel_pos : xy positions in um

    Outputs:
    -------
    cluster_ids : same as input
    is_noise : boolean array, True at index of noise templates

    """
    
    is_noise = np.zeros((templates.shape[0],),dtype='bool')

    print('Checking spread...')
    is_noise += check_template_spread(templates, channel_pos, params)
    print(' Total noise templates: ' + str(np.sum(is_noise)))
    #print(cluster_ids[np.where(is_noise)[0]])

    print('Checking temporal peaks...')
    is_noise += check_template_temporal_peaks(templates, params)
    print(' Total noise templates: ' + str(np.sum(is_noise)))
    #print(cluster_ids[np.where(is_noise)[0]])

    print('Checking spatial peaks...')
    is_noise += check_template_spatial_peaks(templates, channel_pos, params)
    print(' Total noise templates: ' + str(np.sum(is_noise)))
    #print(cluster_ids[np.where(is_noise)[0]])

    return cluster_ids, is_noise[cluster_ids]
    

def check_template_spread(templates, channel_pos, params):

    """
    Checks templates for abnormally large or small channel spread

    Inputs:
    -------
    templates : template for each unit output by Kilosort
    channel_pos : mapping between template channels and actual probe channels

    Outputs:
    -------
    is_noise : boolean array, True at index of noise templates

    Parameters:
    ----------
    """

    is_noise = np.zeros((templates.shape[0],),dtype=bool)

    for i in range(templates.shape[0]):
        MM = np.max(np.abs(templates[i,:,:]),0)
        MM = MM / np.max(MM)
        peak_chan = np.argmax(MM)
        # get channels in this column.
        column_chan = np.where(channel_pos[:,0] == channel_pos[peak_chan,0])
        z_order = np.argsort(channel_pos[column_chan,1])
        z_val = np.sort(channel_pos[column_chan,1])
        MM = MM[column_chan]
        MM = np.squeeze(MM[z_order])   
        
        # calculate column pitch for this column of data        
        z_diff = np.diff(z_val)
        z_diff_unq, counts = np.unique(z_diff, return_counts=True)
        z_pitch = z_diff_unq[np.argmax(counts)]
        
        MMF = gaussian_filter1d(MM, params['smoothed_template_filter_width_um']/z_pitch)

        spread1 = z_pitch * np.sum(MMF > params['smoothed_template_amplitude_threshold'])  # spread of smoothed data
        spread2 = z_pitch * np.sum(MM > params['template_amplitude_threshold'])            # spread of raw data

        # call as noise if:
        # both smoothed and raw data have spread smaller than minimum z spread () or
        # smoothed spread is larger than max
        # characterization of the waveform shape removed for now -- it appears unrelatible across probe types
        
        if (spread1 <= params['mid_spread_threshold_um']):
            is_noise[i] = (spread2 < params['min_spread_threshold_um'])
        elif spread1 > params['max_spread_threshold_um']:
            is_noise[i] = True

    return np.array(is_noise)


def check_template_spatial_peaks(templates, channel_pos, params):

    """
    Checks templates for multiple spatial peaks

    Inputs:
    -------
    templates : template for each unit output by Kilosort
    channel_map : mapping between template channels and actual probe channels

    Outputs:
    -------
    is_noise : boolean array, True at index of noise templates

    Parameters:
    ----------
    """

    nTemplate = templates.shape[0]
    is_noise = np.zeros((nTemplate,),dtype='bool')
    
    # estimate z pitch from all channel_pos
    z_unique = np.unique(channel_pos[:,1])  #returns sorted array
    z_diff = np.diff(z_unique)
    z_diff_unq, counts = np.unique(z_diff, return_counts=True)
    z_pitch = z_diff_unq[np.argmax(counts)]
    
    for index in range(nTemplate):
        is_noise[index] = template_spatial_peaks(templates, channel_pos, z_pitch, params, index)
    # pool = multiprocessing.Pool(np.min([params['multiprocessing_worker_count'],multiprocessing.cpu_count()]))
    # is_noise = pool.map(partial(template_spatial_peaks, templates, channel_map, channel_pos, z_pitch, params), 
    #                     np.arange(templates.shape[0]))

    return np.array(is_noise)


def template_spatial_peaks(templates, channel_pos, z_pitch, params, index):
    
    # JIC notes, 080124
    # altered to use channel postions read from 'channel_positions.npy' file throughout
    # Now only performing interpolation on neightborhood about peak channel, and only
    # for the peak time. Since that is less work, removed running on parallel cpus
    # checks for the presence of alternate peaks of the same sign as the main peak,
    # with amplitude >= main peak * channel_amplitude threshold.

    template = templates[index,:,:]
        
    peak_channel = np.argmax((np.max(template,0) - np.min(template,0)))
    peak_index = np.argmax((np.max(template,1) - np.min(template,1)))
    
    # get template channels within peak_channel_range_um um of peak

    max_dist_um_sq = np.power(params['peak_channel_range_um'],2)
    all_dist_sq = np.power( (channel_pos[:,0] - channel_pos[peak_channel,0]),2) + np.power((channel_pos[:,1] - channel_pos[peak_channel,1]),2)
    chan_inRange = np.argwhere(all_dist_sq < max_dist_um_sq)
    
    ct = np.squeeze(template[:,chan_inRange])
    cp = np.squeeze(channel_pos[chan_inRange,:])


    temp = interpolate_template(ct,cp,z_pitch, peak_index)
    peak_waveform = temp[:,1:6]
    
    pw = peak_waveform.flatten()
    
    si = np.sign(pw[np.argmax(np.abs(pw))])

    peak_locs = []
    
    for x in range(peak_waveform.shape[1]):
        # for row (the interpolated template is reshaped to [nz x nx]
        D = peak_waveform[:,x]
        if np.max(np.abs(D)) >= np.max(np.abs(peak_waveform)) * params['channel_amplitude_thresh']:
            D = D * si
            D = D / np.max(np.abs(D))
            p, peak_prop = find_peaks(D, height = params['peak_height_thresh'], prominence = params['peak_prominence_thresh'])
            if np.any(p):   
                # only count as a 2nd peak if it is 
                if np.argmax(peak_prop['peak_heights']) >= np.max(np.abs(peak_waveform)) * params['channel_amplitude_thresh']:
                    max_ind = np.argmax(peak_prop['peak_heights']) 
                    peak_locs.append(p[max_ind])
              
            
    
    if np.std(peak_locs) >  params['peak_locs_std_thresh']:
        print('Unit id: ' + repr(index))
        print(peak_locs)
        
    return (np.std(peak_locs) > params['peak_locs_std_thresh'])


def check_template_temporal_peaks(templates, params):

    """
    Checks templates for multiple or abnormal temporal peaks

    Inputs:
    -------
    templates : template for each unit output by Kilosort

    Outputs:
    -------
    is_noise : boolean array, True at index of noise templates

    Parameters:
    ----------
    """

    peak_indices = np.argmax((np.max(templates,2) - np.min(templates,2)), 1)

    is_noise = (peak_indices < params['min_temporal_peak_location']) \
               + (peak_indices > params['max_temporal_peak_location'])

    return is_noise



def check_template_shape(template, params):

    """
    Check shape of templates with large spread

    Inputs:
    -------
    template : template for one unit (samples x channels)

    Outputs:
    -------
    is_noise : True if shape is abnormal

    Parameters:
    ----------
    """

    channels_to_use = np.arange(-params['template_shape_channel_range'],
                                params['template_shape_channel_range']+1,
                                4)

    T2 = np.zeros((template.shape[0], channels_to_use.size))
    T2[:] = np.nan

    peak_channel = np.argmax((np.max(template,0) - np.min(template,0)))

    for ii,i in enumerate(channels_to_use):
        try:
            T = template[:,peak_channel+i]
        except IndexError:
            pass
        else:
            T2[:,ii] = T / np.max(np.abs(T))
            

    T3 = T2 - np.tile(T2[:,int(np.floor(channels_to_use.size/2))],
                      (channels_to_use.size,1)
                      ).T
    T4 = np.nanmean(T3,1)
    cwtmatr = cwt(T4, ricker, np.arange(1,template.shape[0],2))
    T5 = cwtmatr[params['wavelet_index'],:]
    wavelet_peak_loc = np.argmax(T5)
    wavelet_peak_height = np.max(T5)

    if wavelet_peak_height > params['min_wavelet_peak_height'] and \
       wavelet_peak_loc > params['min_wavelet_peak_loc'] and \
       wavelet_peak_loc < params['max_wavelet_peak_loc']:
        is_noise = False
    else:
        is_noise = True
    
    return is_noise



def actual_channel_locations(channel_map):
    """
    Physical locations of Neuropixels electrodes, relative to the probe tip
    JIC replaced this inference of position from channel map with 
    direct read of the positions from channel_pos.npy
    
    Inputs:
    -------
    channel_map : mapping between template channels and actual probe channels

    Outputs:
    --------
    locations : (x,y) locations of each electrode (in microns)
    
    """

    max_chan = np.max(channel_map)+1
    actual_channel_locations = np.zeros((max_chan,2))
    xlocations = [16, 48, 0, 32]
    
    for i in range(0, max_chan):
        actual_channel_locations[i,0] = xlocations[i%4]
        actual_channel_locations[i,1] = np.floor(i/2)*20

    return actual_channel_locations[channel_map,:]

def interp_channel_locations(channel_pos, z_pitch):

    """
    Locations of virtual channels after 7x interpolation

    Inputs:
    -------
    channel_pos : channel positions in um

    Outputs:
    --------
    locations : (x,y) locations of each virtual electrode (in microns),
                after 7x interpolation
    
    """

    #max_chan = (np.max(channel_map)+1)*7
    max_chan = channel_pos.shape[0] * 7
    interp_channel_locations = np.zeros((max_chan,2))
    dx = (np.max(channel_pos[:,0]) - np.min(channel_pos[:,0]))/6
    min_x = np.min(channel_pos[:,0])
    xlocations =  min_x + dx*np.asarray(range(7))
    min_z =  np.min(channel_pos[:,1])
    
    for i in range(0, max_chan):
        interp_channel_locations[i,0] = xlocations[i%7]
        interp_channel_locations[i,1] = min_z + np.floor(i/7)*z_pitch/2

    return interp_channel_locations

def interpolate_template(template, channel_pos, z_pitch, peak_index):

    """
    Interpolate template, based on physical channel locations
    JIC altered to read channel positions from output, and limit
    interpolation to a region about peak_channel and peak index

    Inputs:
    -------
    template : template for one unit, wihtin a set distance of the 
               peak_channel (samples x channels)
    channel_positions: positions of channels, in um

    Outputs:
    --------
    template_interp : 3D interpolated template (samples x height x width)
    
    """

    # total_samples = template.shape[0]
    # loc_a = actual_channel_locations(channel_map)
    loc_a = channel_pos             # already correct positions in um
    loc_i = interp_channel_locations(channel_pos, z_pitch)
    
    x_i = np.unique(loc_i[:,0])
    y_i = np.unique(loc_i[:,1])

    interp_temp = np.zeros((len(x_i) * len(y_i)),)
    interp_temp = np.squeeze(griddata(loc_a,np.squeeze(template[peak_index,:]), loc_i, method='cubic', fill_value=0, rescale=False))
   
    return np.reshape(interp_temp, (len(y_i), len(x_i))).astype('float')       

