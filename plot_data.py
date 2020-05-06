import sys
import shutil
import os
import glob
from math import *
import pycrtools as cr
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import datetime

#########################################################################
# Defining functions

def plot_data(file,startsec,endsec,tint,fint,vmin=0.3,vmax=1.2):
    timeresolution=file['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_0/'].attrs['INCREMENT']
    freqaxis=file['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_1/'].attrs['AXIS_VALUES_WORLD']
    startsample=(int(startsec/timeresolution)/tint)*tint
    endsample=(int(endsec/timeresolution)/tint)*tint
    data=file["/SUB_ARRAY_POINTING_000/BEAM_000/STOKES_0"][startsample:endsample,:]
    nfreq=data.shape[1]
    nsamples=endsample-startsample
    data3=np.sum(np.sum(data.reshape(nsamples/tint,tint,nfreq),axis=1).reshape(nsamples/tint,nfreq/fint,fint),axis=2)
    freqaverage=np.average(data3,axis=0)
    plt.imshow((data3/freqaverage).T,origin='lower',aspect='auto',vmin=vmin,vmax=vmax,cmap=plt.get_cmap('hot'),extent=[startsample*timeresolution,endsample*timeresolution,freqaxis[0],freqaxis[-1]])
    return data3

def plot_bf_data(file, plot=True, rficlean=True, cmin=None, cmax=None):

    dynspec = np.load(file)
    
    spectrum = np.median(dynspec,axis=1)
    ds = (dynspec.T/spectrum).T

    ds_median = np.median(dynspec)
    ds_std = np.std(dynspec)
    
    if rficlean:
        # Creating mask
        mask = np.where(spectrum > 5*ds_median)[0]
        not_mask = np.where(spectrum <= 5*ds_median)[0]

        # Fitting spectrum to polynomial
        x = np.arange(len(spectrum))
        p10 = np.poly1d(np.polyfit(x[not_mask], spectrum[not_mask], 10))

        # Updating mask
        mask = np.append(np.where(spectrum > 1.5*p10(x))[0], 
                np.where(spectrum < .5*p10(x))[0])

        for i,s in enumerate(spectrum):
          if i not in mask:
            if len(np.where(dynspec[i] > 2*p10(i))[0]) > dynspec.shape[1]/10:
              #print(i)
              mask = np.concatenate((mask, [i]))

        # Cleaning dynspec
        dynspec_clean = ds
        for i in mask:
            dynspec_clean[i] = np.zeros(ds.shape[1])

        ds_median = np.median(dynspec_clean)
        ds_std = np.std(dynspec_clean)

        for i in range(dynspec.shape[0]):
          if i not in mask:
            filt = np.where(dynspec_clean[i]>1.5*ds_median)[0]
            if len(filt) > dynspec.shape[1]/10:
              mask = np.concatenate((mask, [i]))
              #dynspec_clean[i] = np.zeros(ds.shape[1])
            for j in range(dynspec.shape[1]):
              if dynspec_clean[i,j] > 5*ds_median:
                dynspec_clean[i,j] = ds_median

        not_mask = np.array([i for i in range(len(spectrum)) if i not in mask])

        if plot:
            plt.plot(spectrum, color='teal', label='spectrum')
            plt.fill_between(x,0.5*p10(x),1.5*p10(x), alpha=0.3, 
                    color='turquoise')
            plt.plot(x, p10(x), color='turquoise', 
                    label='allowed spectrum values')
            plt.scatter(mask, spectrum[mask], marker='o', color='orange', 
                    label='mask')

            plt.legend()        
            plt.xlim(0,len(spectrum))
            plt.ylim(1, np.max(spectrum)+1)
            plt.yscale('log')
            plt.show()

    else:
        dynspec_clean = ds
        ds_median = np.median(dynspec_clean)
        ds_std = np.std(dynspec_clean)


    if not cmax: cmax = round((1.4*ds_median), 1)
    if not cmin: cmin = round((0.6*ds_median), 1)
    
    print("Using cmax =", cmax, "cmin =", cmin)

    if plot:
        plt.imshow(dynspec_clean, aspect='auto', origin='lower', 
                interpolation='nearest', cmap='hot', vmin=cmin, vmax=cmax)
        plt.colorbar()
        plt.show()

    # Pulse profile
    dedisp = np.sum(dynspec_clean, axis=0)/(len(spectrum)-len(mask))
    plt.plot(dedisp, color='teal')
    plt.show()
    

#########################################################################

if __name__=='__main__':

    filename=sys.argv[1]
    cmin = sys.argv[2]
    cmax = sys.argv[3]
    #f=h5py.File(filename,'r')
    #d=datetime.datetime.strptime(f.attrs[u'OBSERVATION_START_UTC'][0:19],
    #        '%Y-%m-%dT%H:%M:%S')

    plot_bf_data(filename, rficlean=False, cmin=cmin, cmax=cmax)
