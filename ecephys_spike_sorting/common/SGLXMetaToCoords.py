# -*- coding: utf-8 -*-
"""
Requires python 3
The main() function at the bottom of this file can run from an
interpreter, or, the helper functions can be imported into a
new module or Jupyter notebook.
Standalone program to generate a coordinate file from a SpikeGLX
metadata file for 3A, NP1.0, or NP2 single shank or multishank probes.
The output is set by the outType parameter:
                  0 for text coordinate file; 
                  1 for Kilosort or Kilosort2 channel map file;
                  2 for strings to paste into JRClust .prm file 
                  3 add fields to metadata that were missing prior to SpikeGLX version 032623-phase30
                  4 npy file with (nchan,2) matrix of xy coordinates, e.g. for YASS and related sorters                
@author: Jennifer Colonell, Janelia Research Campus
"""

import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from pathlib import Path
from tkinter import Tk
from tkinter import filedialog
import shutil

# =========================================================
# Parse ini file returning a dictionary whose keys are the metadata
# left-hand-side-tags, and values are string versions of the right-hand-side
# metadata values. We remove any leading '~' characters in the tags to match
# the MATLAB version of readMeta.
#
# The string values are converted to numbers using the "int" and "float"
# fucntions. Note that python 3 has no size limit for integers.
#
def readMeta(metaPath):
    metaDict = {}
    if metaPath.exists():
        # print("meta file present")
        with metaPath.open() as f:
            mdatList = f.read().splitlines()
            # convert the list entries into key value pairs
            for m in mdatList:
                csList = m.split(sep='=')
                if csList[0][0] == '~':
                    currKey = csList[0][1:len(csList[0])]
                else:
                    currKey = csList[0]
                metaDict.update({currKey: csList[1]})
    else:
        print("no meta file")
        
    return(metaDict)
    

# =========================================================
# Return counts of each imec channel type that composes the timepoints
# stored in the binary files.
#
def ChannelCountsIM(meta):
    chanCountList = meta['snsApLfSy'].split(sep=',')
    AP = int(chanCountList[0])
    LF = int(chanCountList[1])
    SY = int(chanCountList[2])
    
    return(AP, LF, SY)
    
# =========================================================
# Return index of SYNC channel in the acquired data
#
def SyncChanIndexIM(meta):
    chanCountList = meta['acqApLfSy'].split(sep=',')
    AP = int(chanCountList[0])
    LF = int(chanCountList[1])
    SY = int(chanCountList[2])
    if SY > 0:
        syInd = AP + LF
    else:    
        syInd = -1
    return(syInd)

# ===========================================================
# Return indicies of the original channels in the acquisition
# In the context of the pipeline, used to build save strings 
# to split out shanks. This list includes the SYNC channel
#
def OriginalChans(meta):
    if meta['snsSaveChanSubset'] == 'all':
        # output = int32, 0 to nSavedChans - 1
        chans = np.arange(0, int(meta['nSavedChans']))
    else:
        # parse the snsSaveChanSubset string
        # split at commas
        chStrList = meta['snsSaveChanSubset'].split(sep=',')
        chans = np.arange(0, 0)  # creates an empty array of int32
        for sL in chStrList:
            currList = sL.split(sep=':')
            if len(currList) > 1:
                # each set of contiguous channels specified by
                # chan1:chan2 inclusive
                newChans = np.arange(int(currList[0]), int(currList[1])+1)
            else:
                newChans = np.arange(int(currList[0]), int(currList[0])+1)
            chans = np.append(chans, newChans)
    return(chans)

# =========================================================
# Return geometry paramters for supported probe types
# These are used to calculate positions from metadata
# that includes only ~snsShankMap
#
def getGeomParams(meta):
# many part numbers have the same geometry parameters ;
# define those sets in lists
# [nShank, shankWidth, shankPitch, even_xOff, odd_xOff, horizPitch, vertPitch, rowsPerShank, elecPerShank]
# offset and pitch values in um
    np1_stag_70um       = [1,  70,   0,   27,   11,  32,  20,  480,  960]
    nhp_lin_70um        = [1,  70,   0,   27,   27,  32,  20,  480,  960]
    nhp_stag_125um_med  = [1, 125,   0,   27,   11,  87,  20, 1368, 2496]
    nhp_stag_125um_long = [1, 125,   0,   27,   11,  87,  20, 2208, 4416]
    nhp_lin_125um_med   = [1, 125,   0,   11,   11, 103,  20, 1368, 2496]
    nhp_lin_125um_long  = [1, 125,   0,   11,   11, 103,  20, 2208, 4416]
    uhd_8col_1bank      = [1,  70,   0,   14,   14,   6,   6,   48,  384]
    uhd_8col_16bank     = [1,  70,   0,   14,   14,   6,   6,  768, 6144]
    np2_ss              = [1,  70,   0,   27,   27,  32,  15,  640, 1280]
    np2_4s              = [4,  70, 250,   27,   27,  32,  15,  640, 1280]
    NP1120              = [1,  70,   0, 6.75, 6.75, 4.5, 4.5,  192,  384]
    NP1121              = [1,  70,   0, 6.25, 6.25,   3,   3,  384,  384]
    NP1122              = [1,  70,   0, 6.75, 6.75, 4.5, 4.5,   24,  384]
    NP1123              = [1,  70,   0,10.25,10.25,   3,   3,   32,  384]
    NP1300              = [1,  70,   0,   11,   11,  48,  20,  480,  960]
    NP1200              = [1,  70,   0,   27,   11,  32,  20,   64,  128]
    NXT3000             = [1,  70,   0,   53,   53,   0,  15,  128,  128]
    
    M = dict([
    ('3A',np1_stag_70um),
    ('PRB_1_4_0480_1',np1_stag_70um),
    ('PRB_1_4_0480_1_C',np1_stag_70um),
    ('NP1010',np1_stag_70um), 
    ('NP1011',np1_stag_70um),
    ('NP1012',np1_stag_70um),
    ('NP1013',np1_stag_70um),
    
    ('NP1015',nhp_lin_70um),
    ('NP1015',nhp_lin_70um),
    ('NP1016',nhp_lin_70um),
    ('NP1017',nhp_lin_70um),
       
    ('NP1020',nhp_stag_125um_med),
    ('NP1021',nhp_stag_125um_med),
    ('NP1030',nhp_stag_125um_long),
    ('NP1031',nhp_stag_125um_long),
    
    ('NP1022',nhp_lin_125um_med),
    ('NP1032',nhp_lin_125um_long),
    
    ('NP1100',uhd_8col_1bank),
    ('NP1110',uhd_8col_16bank),
    
    ('PRB2_1_4_0480_1',np2_ss),
    ('PRB2_1_2_0640_0',np2_ss),
    ('NP2000',np2_ss),
    ('NP2003',np2_ss),
    ('NP2004',np2_ss),
    
    ('PRB2_4_2_0640_0',np2_4s),
    ('PRB2_4_4_0480_1',np2_4s),
    ('NP2010',np2_4s),
    ('NP2013',np2_4s),
    ('NP2014',np2_4s),
    
    ('NP1120',NP1120),
    ('NP1121',NP1121),
    ('NP1122',NP1122),
    ('NP1123',NP1123),
    ('NP1300',NP1300),
    
    ('NP1200',NP1200),
    ('NXT3000',NXT3000)
    ])
    
    # get probe part number; if absent, this is a 3A
    if 'imDatPrb_pn' in meta:   
        pn = meta['imDatPrb_pn']
    else:
        pn = '3A';


    if pn in M:
        geomList = M[pn]
    else:
        print('unsupported probe part number\n')
        geomList = []
    
    return geomList


# =========================================================
# Return full MUX table string to append to metadata
#
def getMuxTable(meta):
# Read probe part number from meta
# Return full MUX table string to append to metadata
# As of 032923, there are 4 mux tables

    np1=r'~muxTbl=(32,12)(0 1 24 25 48 49 72 73 96 97 120 121 144 145 168 169 192 193 216 217 240 241 264 265 288 289 312 313 336 337 360 361)(2 3 26 27 50 51 74 75 98 99 122 123 146 147 170 171 194 195 218 219 242 243 266 267 290 291 314 315 338 339 362 363)(4 5 28 29 52 53 76 77 100 101 124 125 148 149 172 173 196 197 220 221 244 245 268 269 292 293 316 317 340 341 364 365)(6 7 30 31 54 55 78 79 102 103 126 127 150 151 174 175 198 199 222 223 246 247 270 271 294 295 318 319 342 343 366 367)(8 9 32 33 56 57 80 81 104 105 128 129 152 153 176 177 200 201 224 225 248 249 272 273 296 297 320 321 344 345 368 369)(10 11 34 35 58 59 82 83 106 107 130 131 154 155 178 179 202 203 226 227 250 251 274 275 298 299 322 323 346 347 370 371)(12 13 36 37 60 61 84 85 108 109 132 133 156 157 180 181 204 205 228 229 252 253 276 277 300 301 324 325 348 349 372 373)(14 15 38 39 62 63 86 87 110 111 134 135 158 159 182 183 206 207 230 231 254 255 278 279 302 303 326 327 350 351 374 375)(16 17 40 41 64 65 88 89 112 113 136 137 160 161 184 185 208 209 232 233 256 257 280 281 304 305 328 329 352 353 376 377)(18 19 42 43 66 67 90 91 114 115 138 139 162 163 186 187 210 211 234 235 258 259 282 283 306 307 330 331 354 355 378 379)(20 21 44 45 68 69 92 93 116 117 140 141 164 165 188 189 212 213 236 237 260 261 284 285 308 309 332 333 356 357 380 381)(22 23 46 47 70 71 94 95 118 119 142 143 166 167 190 191 214 215 238 239 262 263 286 287 310 311 334 335 358 359 382 383)'
    np2=r'~muxTbl=(24,16)(0 1 32 33 64 65 96 97 128 129 160 161 192 193 224 225 256 257 288 289 320 321 352 353)(2 3 34 35 66 67 98 99 130 131 162 163 194 195 226 227 258 259 290 291 322 323 354 355)(4 5 36 37 68 69 100 101 132 133 164 165 196 197 228 229 260 261 292 293 324 325 356 357)(6 7 38 39 70 71 102 103 134 135 166 167 198 199 230 231 262 263 294 295 326 327 358 359)(8 9 40 41 72 73 104 105 136 137 168 169 200 201 232 233 264 265 296 297 328 329 360 361)(10 11 42 43 74 75 106 107 138 139 170 171 202 203 234 235 266 267 298 299 330 331 362 363)(12 13 44 45 76 77 108 109 140 141 172 173 204 205 236 237 268 269 300 301 332 333 364 365)(14 15 46 47 78 79 110 111 142 143 174 175 206 207 238 239 270 271 302 303 334 335 366 367)(16 17 48 49 80 81 112 113 144 145 176 177 208 209 240 241 272 273 304 305 336 337 368 369)(18 19 50 51 82 83 114 115 146 147 178 179 210 211 242 243 274 275 306 307 338 339 370 371)(20 21 52 53 84 85 116 117 148 149 180 181 212 213 244 245 276 277 308 309 340 341 372 373)(22 23 54 55 86 87 118 119 150 151 182 183 214 215 246 247 278 279 310 311 342 343 374 375)(24 25 56 57 88 89 120 121 152 153 184 185 216 217 248 249 280 281 312 313 344 345 376 377)(26 27 58 59 90 91 122 123 154 155 186 187 218 219 250 251 282 283 314 315 346 347 378 379)(28 29 60 61 92 93 124 125 156 157 188 189 220 221 252 253 284 285 316 317 348 349 380 381)(30 31 62 63 94 95 126 127 158 159 190 191 222 223 254 255 286 287 318 319 350 351 382 383)'
    np1100=r'~muxTbl=(32,12)(0 1 24 25 48 49 72 73 96 97 120 121 144 145 168 169 192 193 216 217 240 241 264 265 288 289 312 313 336 337 360 361)(2 3 26 27 50 51 74 75 98 99 122 123 146 147 170 171 194 195 218 219 242 243 266 267 290 291 314 315 338 339 362 363)(4 5 28 29 52 53 76 77 100 101 124 125 148 149 172 173 196 197 220 221 244 245 268 269 292 293 316 317 340 341 364 365)(6 7 30 31 54 55 78 79 102 103 126 127 150 151 174 175 198 199 222 223 246 247 270 271 294 295 318 319 342 343 366 367)(8 9 32 33 56 57 80 81 104 105 128 129 152 153 176 177 200 201 224 225 248 249 272 273 296 297 320 321 344 345 368 369)(10 11 34 35 58 59 82 83 106 107 130 131 154 155 178 179 202 203 226 227 250 251 274 275 298 299 322 323 346 347 370 371)(12 13 36 37 60 61 84 85 108 109 132 133 156 157 180 181 204 205 228 229 252 253 276 277 300 301 324 325 348 349 372 373)(14 15 38 39 62 63 86 87 110 111 134 135 158 159 182 183 206 207 230 231 254 255 278 279 302 303 326 327 350 351 374 375)(16 17 40 41 64 65 88 89 112 113 136 137 160 161 184 185 208 209 232 233 256 257 280 281 304 305 328 329 352 353 376 377)(18 19 42 43 66 67 90 91 114 115 138 139 162 163 186 187 210 211 234 235 258 259 282 283 306 307 330 331 354 355 378 379)(20 21 44 45 68 69 92 93 116 117 140 141 164 165 188 189 212 213 236 237 260 261 284 285 308 309 332 333 356 357 380 381)(22 23 46 47 70 71 94 95 118 119 142 143 166 167 190 191 214 215 238 239 262 263 286 287 310 311 334 335 358 359 382 383)'
    np128ch=r'~muxTbl=(12,12)(84 11 85 5 74 10 56 112 46 121 39 127)(100 26 110 33 69 24 63 109 45 93 25 99)(87 0 82 6 71 15 53 117 43 122 42 116)(102 28 81 34 70 18 60 103 17 94 27 101)(73 1 86 7 68 16 50 106 40 123 128 128)(105 29 75 35 67 12 54 89 20 95 128 128)(76 2 83 8 65 13 47 118 49 124 128 128)(108 30 78 36 64 14 51 90 23 96 128 128)(79 3 80 9 62 114 44 119 52 125 128 128)(104 31 72 37 61 113 57 91 19 97 128 128)(88 4 77 21 59 111 41 120 55 126 128 128)(107 32 66 38 58 115 48 92 22 98 128 128)'
       
    M = dict([
    ('3A',np1),
    ('PRB_1_4_0480_1',np1),
    ('PRB_1_4_0480_1_C',np1),
    ('NP1010',np1), 
    ('NP1011',np1),
    ('NP1012',np1),
    ('NP1013',np1),
    
    ('NP1015',np1),
    ('NP1015',np1),
    ('NP1016',np1),
    ('NP1017',np1),
       
    ('NP1020',np1),
    ('NP1021',np1),
    ('NP1030',np1),
    ('NP1031',np1),
    
    ('NP1022',np1),
    ('NP1032',np1),
    
    ('NP1100',np1),
    ('NP1110',np1100),
    
    ('PRB2_1_4_0480_1',np2),
    ('PRB2_1_2_0640_0',np2),
    ('NP2000',np2),
    ('NP2003',np2),
    ('NP2004',np2),
    
    ('PRB2_4_2_0640_0',np2),
    ('PRB2_4_4_0480_1',np2),
    ('NP2010',np2),
    ('NP2013',np2),
    ('NP2014',np2),
    
    ('NP1120',np1),
    ('NP1121',np1),
    ('NP1122',np1),
    ('NP1123',np1),
    ('NP1300',np1),
    
    ('NP1200',np128ch),
    ('NXT3000',np128ch)
    ])
    
    # get probe part number; if absent, this is a 3A
    if 'imDatPrb_pn' in meta:   
        pn = meta['imDatPrb_pn']
    else:
        pn = '3A';


    if pn in M:
        muxTableStr = M[pn] + '\n'
    else:
        print('unsupported probe part number\n')
        muxTableStr = []
    
    return muxTableStr


# =========================================================
# Parse imro table to extract 'new' metadata items
# (SpikeGLX 032623-phase 30 and later)
# Returns strings for:
#     imChan0apGain
#     imChan0lfGain
#     imAnyChanFullBand
# 
#
def imroMetaItems(meta):
        
    # read in the imro table 
    imroTbl = meta['imroTbl'].split(sep=')')
    
    # there is an entry in the map for each saved channel
    # number of entries in map -- subtract 1 for header, one for trailing ')'
    nEntry = len(imroTbl) - 2

    # read header to get what type of imro this is
    headEntry = imroTbl[0]
    headEntry = headEntry[1:len(headEntry)]   #remove leading '('
    headList = headEntry.split(sep=',')
    currType = headList[0]
    
    if int(currType) > 50000:
        # this is a 3A probe
        currType = '0'

    if currType == '24' or currType == '21':
        imChan0apGainStr = 'imChan0apGain=80'
        imChan0lfGainStr = 'imChan0lfGain=80'
        imAnyChanFullBandStr = 'imAnyChanFullBand=true'
        
    elif currType == '0':
        # parse first entry of the table for lf and apgain
        currEntry = imroTbl[1]
        currEntry = currEntry[1:len(currEntry)]   #remove leading '('
        currList = currEntry.split(sep=' ')
        imChan0apGainStr = 'imChan0apGain=' + currList[3]
        imChan0lfGainStr = 'imChan0lfGain=' + currList[4]
                   
        imAnyChanFullBandStr = 'imAnyChanFullBand=false'
        if len(currList) == 6:  # indicates this is not a 3A imro table
        #check for any full band
            for i in range(nEntry):
                # check for any channels where the AP filter was not used
                currEntry = imroTbl[i+1]
                currEntry = currEntry[1:len(currEntry)]
                currList = currEntry.split(sep=' ')
                if currList[5] == '0':
                    imAnyChanFullBandStr = 'imAnyChanFullBand=true'
                    break

    elif currType == '1110':
        # ap and lf gain and filter option are in the header 
        imChan0apGainStr = 'imChan0apGain=' + headList[3]
        imChan0lfGainStr = 'imChan0lfGain=' + headList[4]
        imAnyChanFullBandStr = 'imAnyChanFullBand=false'
        if headList[5] == 0:
            imAnyChanFullBandStr = 'imAnyChanFullBand=true'
   
    
    # add newline characters to each
    imChan0apGainStr = imChan0apGainStr + '\n'
    imChan0lfGainStr = imChan0lfGainStr + '\n'
    imAnyChanFullBandStr = imAnyChanFullBandStr + '\n'
    
    return imChan0apGainStr, imChan0lfGainStr, imAnyChanFullBandStr


# =========================================================
# Parse snsGeomMap for XY coordinates
#
def geomMapToGeom(meta):

    # read in the shank map   
    geomMap = meta['snsGeomMap'].split(sep=')')
    
    # there is an entry in the map for each saved channel
    # number of entries in map -- subtract 1 for header, one for trailing ')'
    nEntry = len(geomMap) - 2

    shankInd = np.zeros((nEntry,))
    xCoord = np.zeros((nEntry,))
    yCoord = np.zeros((nEntry,))
    connected = np.zeros((nEntry,))
    
    for i in range(nEntry):
        # get parameter list from this entry, skipping first header entry
        currEntry = geomMap[i+1]
        currEntry = currEntry[1:len(currEntry)]
        currList = currEntry.split(sep=':')
        shankInd[i] = int(currList[0])
        xCoord[i] = float(currList[1])
        yCoord[i] = float(currList[2])
        connected[i] = int(currList[3])

    # parse header for number of shanks
    currList = geomMap[0].split(',')
    nShank = int(currList[1]);
    shankPitch = float(currList[2]);
    shankWidth = float(currList[3]);
    
    return nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected


# =========================================================
# Build snsGeomMap from xy coordinates
#
def snsGeom(meta, shankInd, xCoord, yCoord, use):
    # header
    # get probe part number; if absent, this is a 3A
    if 'imDatPrb_pn' in meta:   
        pn = meta['imDatPrb_pn']
    else:
        pn = '3A';
        
    geomList = getGeomParams(meta)
    # geomList =
    # [nShank, shankWidth, shankPitch, even_xOff, odd_xOff, horizPitch, vertPitch, rowsPerShank, elecPerShank]
    
    snsGeomStr = '~snsGeomMap=(' + pn
    snsGeomStr = snsGeomStr + ',{:d},{:g},{:g})'.format(geomList[0], geomList[2], geomList[1])
    nEntry = shankInd.shape[0]
    
    for i in range(0,nEntry):
        snsGeomStr = snsGeomStr + '({:g}:{:g}:{:g}:{:g})'.format(shankInd[i], xCoord[i], yCoord[i], use[i])
    
    snsGeomStr = snsGeomStr + '\n'
    
    return snsGeomStr


# =========================================================
# Get XY coordinates from snsShankMap plus hard coded geom values
#
def shankMapToGeom(meta):
    # get number of saved AP channels (some early metadata files have a
    # SYNC entry in the snsChanMap)
    AP, LF, SY = ChannelCountsIM(meta)

    shankMap = meta['snsShankMap'].split(sep=')')
    
    shankInd = np.zeros((AP,))
    colInd = np.zeros((AP,))
    rowInd = np.zeros((AP,))
    connected = np.zeros((AP,))
    xCoord = np.zeros((AP,))
    yCoord = np.zeros((AP,))
    
    for i in range(AP):
        # get parameter list from this entry, skipping first header entry
        currEntry = shankMap[i+1]
        currEntry = currEntry[1:len(currEntry)]
        currList = currEntry.split(sep=':')
        shankInd[i] = int(currList[0])
        colInd[i] = int(currList[1])
        rowInd[i] = int(currList[2])
        connected[i] = int(currList[3])
   
    geomList = getGeomParams(meta);
    # geomList = 
    # [nShank, shankWidth, shankPitch, even_xOff, odd_xOff, horizPitch, vertPitch, rowsPerShank, elecPerShank]
    
    oddRows = np.bool_(rowInd%2);
    evenRows = ~oddRows;
    xCoord = colInd*float(geomList[5]);
    xCoord[evenRows] = xCoord[evenRows] + geomList[3]
    xCoord[oddRows] = xCoord[oddRows] + geomList[4]
    yCoord = rowInd*float(geomList[6])
    
    nShank = geomList[0]
    shankWidth = geomList[1]
    shankPitch = geomList[2]
    
    return nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected


# =========================================================
# Plot x z positions of all electrodes and saved channels
#
def plotSaved(xCoord, yCoord, shankInd, meta):
    
    geomList = getGeomParams(meta)
    # geomList = 
    # [nShank, shankWidth, shankPitch, even_xOff, odd_xOff, horizPitch, vertPitch, rowsPerShank, elecPerShank]
    
    # calculate positions on one shank
    nCol = geomList[8]/geomList[7]
    rowInd = np.arange(geomList[8])
    rowInd = np.floor(rowInd/nCol)
    oddRows = np.bool_(rowInd%2);
    evenRows = ~oddRows;
    
    colInd = np.arange(geomList[8])
    colInd = (colInd % nCol)
    
    xall = colInd*geomList[5]
    xall[evenRows] = xall[evenRows] + geomList[3]
    xall[oddRows] = xall[oddRows] + geomList[4]
    
    yall = rowInd*geomList[6];
    
    
    fig = plt.figure(figsize=(2,12))
    
    shankSep = geomList[2]
    
    # loop over shanks
    for sI in range(geomList[0]):
        
        # plot all positions   
        marker_style = dict(c='w', edgecolor = 'k', linestyle='None', marker='s', s=5)                 
        plt.scatter(shankSep*sI + xall, yall, **marker_style)
        
        # plot selected positions
        currInd = np.argwhere(shankInd == sI)
        marker_style = dict(c='b', edgecolor = 'g', linestyle='None', marker='s', s=15) 
        plt.scatter(shankSep*sI + xCoord[currInd], yCoord[currInd], **marker_style)
   
   # after looping over all shanks, show plot 
    plt.show()

    return


# =========================================================
# CoordsTo... functions to write output in different 
# formats
#
def CoordsToText(meta, chans, xCoord, yCoord, connected, shankInd, shankSep, baseName, savePath, buildPath ):

    if buildPath:
        newName = baseName + '_siteCoords.txt'
        saveFullPath = Path(savePath / newName)
    else:
        saveFullPath = savePath

    # Note that the channel index written is the index of that channel in the saved file

    with open(saveFullPath, 'w') as outFile:
        for i in range(0,chans.size):
            currX = shankInd[i]*shankSep + xCoord[i]
            currLine = '{:d}\t{:g}\t{:g}\t{:g}\n'.format(i, currX, yCoord[i], shankInd[i])
            outFile.write(currLine)


def CoordsToNPY(meta, chans, xCoord, yCoord, connected, shankInd, shankSep, baseName, savePath, buildPath ):

    if buildPath:
        newName = baseName + '_siteCoords.npy'
        saveFullPath = Path(savePath / newName)
    else:
        saveFullPath = savePath

    # write an npy file of nChanx2 
    nchan = xCoord.shape[0]
    geom = np.zeros((nchan,2))
    geom[:,0] = xCoord + shankInd*shankSep
    geom[:,1] = yCoord
    
    np.save(saveFullPath, geom)

            
def CoordsToJRCString(meta, chans, xCoord, yCoord, connected, shankInd, shankSep, baseName, savePath, buildPath ):
    
    if buildPath:
        newName = baseName +'_forJRCprm.txt'
        saveFullPath = Path(savePath / newName)
    else:
        saveFullPath = savePath
    
    # siteMap, equivalent of chanMap in KS, is the order of channels in the saved file
    # rather than original channel indicies.
    nChan = chans.size
    siteMap = np.arange(0,nChan, dtype = 'int')
    siteMap = siteMap + 1   # convert to 1-based for MATLAB
    
    shankInd = shankInd + 1 # conver to 1-based for MATLAB
 
    shankStr = 'shankMap = ['
    coordStr = 'siteLoc = ['
    siteMapStr = 'siteMap = ['
    
    xCoord = shankInd*shankSep + xCoord
    
    for i in range(0,chans.size-1):
        shankStr = shankStr + '{:g},'.format(shankInd[i])   # convert to 1-based for MATLAB
        coordStr = coordStr + '{:g},{:g};'.format(xCoord[i], yCoord[i])
        siteMapStr = siteMapStr + '{:d},'.format(siteMap[i])    # convert to 1-based for MATLAB
    
    # final entries
    shankStr = shankStr + '{:g}];\n'.format(shankInd[nChan-1])
    coordStr = coordStr + '{:g},{:g}];\n'.format(xCoord[nChan-1], yCoord[nChan-1])
    siteMapStr = siteMapStr + '{:d}];\n'.format(siteMap[nChan-1])
    
    with open(saveFullPath, 'w') as outFile:
        outFile.write(shankStr) 
        outFile.write(coordStr)
        outFile.write(siteMapStr)


def CoordsToKSChanMap(meta, chans, xCoord, yCoord, connected, shankInd, shankSep, baseName, savePath, buildPath ): 
    
    if buildPath:
        newName = baseName +'_kilosortChanMap.mat'
        saveFullPath = Path(savePath / newName)
    else:
        saveFullPath = savePath
    
    nChan = chans.size
    # channel map is the order of channels in the file, rather than the 
    # original indicies of the channels
    chanMap0ind = np.arange(0,nChan,dtype='float64')
    chanMap0ind = chanMap0ind.reshape((nChan,1))
    chanMap = chanMap0ind + 1
    
    connected = (connected==1)
    connected = connected.reshape((nChan,1))
    
    xCoord = shankInd*shankSep + xCoord
    xCoord = xCoord.reshape((nChan,1))
    yCoord = yCoord.reshape((nChan,1))
    
    kcoords = shankInd + 1
    kcoords = kcoords.reshape((nChan,1))
    kcoords = kcoords.astype('float64')
    
    name = baseName
    
    mdict = {
            'chanMap':chanMap,
            'chanMap0ind':chanMap0ind,
            'connected':connected,
            'name':name,
            'xcoords':xCoord,
            'ycoords':yCoord,
            'kcoords':kcoords,
            }
    scipy.io.savemat(saveFullPath, mdict)

    
def CoordsToGeomMap(meta, chans, xCoord, yCoord, connected, shankInd, shankSep, baseName, savePath, buildPath ):

    if buildPath:
        newName = baseName + '_orig.meta'
        copyFullPath = Path(savePath / newName)
    else:
        print('Can only make new metadata in same directory as original.')
        return
    
    origPath = Path(savePath / (baseName + '.meta'))
    shutil.move(origPath, copyFullPath)
    shutil.copy(copyFullPath,origPath)
    
    # check 'new' fields; add if not present, add everything (most common case)
    if 'imChan0apGain' not in meta:        
        imChan0apGainStr, imChan0lfGainStr, imAnyChanFullBandStr = imroMetaItems(meta)
        muxTableStr = getMuxTable(meta)        
        snsGeomStr = snsGeom(meta, shankInd, xCoord, yCoord, connected)
    
        with open(origPath, 'a') as outFile:
            outFile.write(imChan0apGainStr)
            outFile.write(imChan0lfGainStr)
            outFile.write(imAnyChanFullBandStr)
            outFile.write(muxTableStr)
            outFile.write(snsGeomStr)
            


# =========================================================
# Given a path to a SpikeGLX metadata file, write out coordinates
# in formats for analysis software to consume
# Input params:
#   metaFullPath: full path, including the file name
#   outType:  format for the output
#   badChan:  channels other than reference channels to exclude
#   destFullPath: 
#  
def MetaToCoords(metaFullPath, outType, badChan= np.zeros((0), dtype = 'int'), destFullPath = '', showPlot=False):
       
    # Read in metadata; returns a dictionary with string for values
    meta = readMeta(metaFullPath)
    
    # Get coordinates for saved channels from snsGeomMap, if present,
    # otherwise from snsShankMap
    if 'snsGeomMap' in meta:
        [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(meta)
    else:   
        [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = shankMapToGeom(meta)
    
    if showPlot:
        plotSaved(xCoord, yCoord, shankInd, meta)
        
    [AP,LF,SY] = ChannelCountsIM(meta)    
    chans = np.arange(AP)
    
    # Channels identified as noisy by kilosort helper indexed
    # according to position in the file
    # since these can include the SYNC channel, remove any from
    # list that are outside the range of AP channels
    badChan = badChan[badChan < AP]
    connected[badChan] = 0    

    baseName = metaFullPath.stem
    
    if outType >= 0:
        if len(destFullPath) == 0:
            savePath = metaFullPath.parent
            buildPath = True
        else:
            buildPath = False
            savePath = destFullPath
        outputSwitch = {
                0: CoordsToText,
                1: CoordsToKSChanMap,
                2: CoordsToJRCString,
                3: CoordsToGeomMap,
                4: CoordsToNPY
        }
        
        writeFunc = outputSwitch.get(outType)
        writeFunc(meta, chans, xCoord, yCoord, connected, shankInd, shankPitch, baseName, savePath, buildPath )
    
    return xCoord, yCoord, shankInd, connected


# =========================================================    
# Sample calling program to get a metadata file from the user,
# output a file set by outType
#   0 = tab delimited text file of coordinates in um, index, x, y, shank index
#   1 = KS2 chan map .mat file
#   2 = strings of channel map, shank index, and x,y pairs for JRClust
#   3 = create a new metadata file with snsGeomMap appended (for converting old metadata to new)
#   4 = npy file with (nchan,2) matrix of xy coordinates, e.g. for YASS and related sorters
# file is saved to the same path as the metadata file.
#
def main():
    
    outType = 1
    
    # Get file from user
    root = Tk()         # create the Tkinter widget
    root.withdraw()     # hide the Tkinter root window

    # Windows specific; forces the window to appear in front
    root.attributes("-topmost", True)

    metaFullPath = Path(filedialog.askopenfilename(title="Select meta file"))
    root.destroy()      # destroy the Tkinter widget
   
    MetaToCoords( metaFullPath=metaFullPath, outType=outType, showPlot=True  )

        

if __name__ == "__main__":
    main()