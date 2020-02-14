import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import datetime

filename=sys.argv[1]
f=h5py.File(filename,'r')
d=datetime.datetime.strptime(f.attrs[u'OBSERVATION_START_UTC'][0:19],'%Y-%m-%dT%H:%M:%S')
timeresolution=f['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_0/'].attrs['INCREMENT']
freqaxis=f['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_1/'].attrs['AXIS_VALUES_WORLD']

def dedisperse(data, dm, dt=0.00032768, freq=(100., 200.), freq_ref=150.0):
    data = data.copy()
    
    nfreq, ntime = data.shape[0], data.shape[1]

    freqs = np.linspace(freq[0], freq[-1], nfreq)

    if freq_ref is None:
        freq_ref = freqs.max()

    tdelay = 4.148e3*dm*(freqs**-2 - freq_ref**-2)
    ntime = len(data[0])

    maxind_arr = []

    for ii, f in enumerate(freqs):
        data[ii] = np.roll(data[ii], -np.int(tdelay[ii]/dt))

    return data

def plot_data(file,startsec,endsec,tint,fint,vmin=-3,vmax=3, rficlean=True, timeresolution=0.00032768):
    try:
        timeresolution=file['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_0/'].attrs['INCREMENT']
        freqaxis=file['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_1/'].attrs['AXIS_VALUES_WORLD']
        startsample=(int(startsec/timeresolution)/tint)*tint
        endsample=(int(endsec/timeresolution)/tint)*tint
        data=file["/SUB_ARRAY_POINTING_000/BEAM_000/STOKES_0"][startsample:endsample,:]
    except:
        data = file
        startsample=(int(startsec/timeresolution)/tint)*tint
        endsample=(int(endsec/timeresolution)/tint)*tint

    if rficlean:
        print('RFI Cleaning')
        sig = np.std(data, axis=0)
        ind_mask = np.where(sig>np.median(sig)*3)[0]
        data[:, ind_mask] = 0.0

    nfreq=data.shape[1]
    nsamples=endsample-startsample
    data = data[:nsamples//tint*tint, :nfreq//fint*fint]
    data3=np.sum(np.sum(data.reshape(nsamples/tint,tint,nfreq),axis=1).reshape(nsamples/tint,nfreq/fint,fint),axis=2)
    freqaverage=np.average(data3,axis=0)
    data3 -= np.median(data3, axis=0, keepdims=True)
    data3 /= np.std(data3, axis=0, keepdims=True)
    data3[np.isnan(data3)] = 0.0
    data3 = data3.T
    plt.imshow(data3, origin='lower',aspect='auto',vmin=vmin,vmax=vmax,cmap=plt.get_cmap('hot'),
            extent=[startsample*timeresolution,endsample*timeresolution,freqaxis[0]*1e-6,freqaxis[-1]*1e-6])
    plt.ylabel('Freq [Mhz]')
    plt.xlabel('Time [s]')
    plt.show()
    return data3

def plot_dedispersed_timestream(filename, startsec,endsec,tint,fint, dm=26.76, rficlean=True):

    try:
        file=h5py.File(filename,'r')
        d=datetime.datetime.strptime(file.attrs[u'OBSERVATION_START_UTC'][0:19],'%Y-%m-%dT%H:%M:%S')
        timeresolution=file['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_0/'].attrs['INCREMENT']
        freqaxis=file['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_1/'].attrs['AXIS_VALUES_WORLD']
        startsample=(int(startsec/timeresolution)/tint)*tint
        endsample=(int(endsec/timeresolution)/tint)*tint
        data=file["/SUB_ARRAY_POINTING_000/BEAM_000/STOKES_0"][startsample:endsample,:]
    except:
        data = file
        startsample=(int(startsec/timeresolution)/tint)*tint
        endsample=(int(endsec/timeresolution)/tint)*tint

    if rficlean:
        print('RFI Cleaning')
        sig = np.std(data, axis=0)
        ind_mask = np.where(sig>np.median(sig)*3)[0]
        data[:, ind_mask] = 0.0

    nfreq=data.shape[1]
    nsamples=endsample-startsample
    data = data[:nsamples//tint*tint, :nfreq//fint*fint]
    data3=np.sum(np.sum(data.reshape(nsamples/tint,tint,nfreq),axis=1).reshape(nsamples/tint,nfreq/fint,fint),axis=2)
    nsamples = data3.shape[0]
    freqaverage=np.average(data3,axis=0)
    data3 -= np.median(data3, axis=0, keepdims=True)
    data3 /= np.std(data3, axis=0, keepdims=True)
    data3[np.isnan(data3)] = 0.0
    timeresolution *= tint 
    times_arr = np.linspace(startsec,endsec,nsamples)

    data_dedisp = dedisperse(data3.T, dm, dt=timeresolution, freq=(freqaxis[0], freqaxis[-1]), freq_ref=150.0)
    plt.plot(times_arr, data_dedisp.mean(0))
    plt.xlabel('Time [s]', fontsize=20)
    plt.show()

    return data3, data_dedisp, timeresolution



print "determine starttime"
print "1574763544-int(d.strftime('%s'))"
print "plot_data(f,1140,1148,4,8,0.2,1.6)"

#plot_data(f,100,105,4,8,3,-2)




