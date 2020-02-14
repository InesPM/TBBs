import h5py
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import argparse
import dedispersion as dd

#-------------------------------------------------------------------------------
# Defining functions
#-------------------------------------------------------------------------------

def data(file, DM=26.8, nch=16, t_int=6):
    """
    Function extracting dynamic spectrum from input file
    Inputs:
    file  - *.h5 input file name
    DM    - dispersion measure of the source
    nch   - number of frequency channels
    t_int - time channels
    """
    f_bf=h5py.File(file,'r')

    st=list(f_bf.keys())[0]
    # We need to make sure the subbands are read in order
    subbands=sorted(list(f_bf[st].keys()))
    nsb=len(subbands)
    datalen=f_bf[st][subbands[0]].shape[0]
    datalen=datalen-datalen%(2*nch*t_int)

    complexdynspec=np.zeros(dtype=f_bf[st][subbands[0]].dtype,shape=(nsb,datalen))
    subband_frequencies=np.zeros(dtype=float,shape=(nsb))

    for num,sb in enumerate(subbands):
        complexdynspec[num]=f_bf[st][sb][0:datalen]
        subband_frequencies[num]=f_bf[st][sb].attrs[u'CENTRAL_FREQUENCY']

    # Same resolution as BF data : nch=16
    data01=complexdynspec[:,:].T
    s=data01.shape; fftdata01=np.fft.fft(data01.reshape(s[0]//nch,nch,s[1]).swapaxes(1,2),axis=2).reshape(s[0]//nch,nch*s[1])

    # Dynamic spectrum
    assert(f_bf[st][subbands[0]].attrs['BANDWIDTH_UNIT']==b'Hz'), "Bandwidth not in Hz"
    subband_bandwidth=f_bf[st][subbands[0]].attrs['BANDWIDTH'] #Hz
    assert(f_bf[st][subbands[0]].attrs['TIME_RESOLUTION_UNIT']==b's'), "Time resolution not in s"
    subband_timeresolution=f_bf[st][subbands[0]].attrs['TIME_RESOLUTION'] #seconds

    channel_timeresolution=subband_timeresolution*nch
    channel_bandwidth=subband_bandwidth/nch
    channel_frequencies=np.zeros(dtype=float,shape=(nsb,nch))
    channel_delays=np.zeros(dtype=int,shape=(nsb,nch))

    for num,sbfreq in enumerate(subband_frequencies):
        channel_frequencies[num]=(np.arange(nch)-nch/2)*channel_bandwidth+sbfreq
        channel_delays[num]=dd.freq_to_delay(DM,channel_frequencies[num],channel_timeresolution)
    channel_frequencies=channel_frequencies.reshape((nsb*nch))
    channel_delays=channel_delays.reshape((nsb*nch))

    dynspec1=np.abs(fftdata01)
    if False:
        for ch,delay in enumerate(channel_delays):
            if int(delay)>=1:
                dynspec1[0:-int(delay),ch]=dynspec1[int(delay):,ch]
                dynspec1[-int(delay):,ch]=0 #dynspec[0:int(delay),ch]
    #print("Dynspec shape : ", dynspec.shape)

    # Spectrum
    specstd=np.std(dynspec1,axis=0)
    spect_raw=np.sum(dynspec1,axis=0)
    noisyness=(spect_raw/np.median(spect_raw))/(specstd/np.median(specstd))
    flagfreq=(noisyness<0.95)|(noisyness>1.05)

    # Start time
    starttime=[f_bf[st][sb].attrs[u'TIME']%100+f_bf[st][sb].attrs[u'SLICE_NUMBER']*subband_timeresolution for sb in sorted(subbands)]

    # Downsampling data
    f_int=nch//2 # Best visualized when nch/f_int=2
    dynspec2=np.sum(dynspec1.reshape(dynspec1.shape[0]//t_int,t_int,dynspec1.shape[1])/t_int,axis=1)
    dynspec3=np.sum(dynspec2.reshape(dynspec2.shape[0],dynspec2.shape[1]//f_int,f_int)/f_int,axis=2)
    spectrum=np.average(dynspec3,axis=0)

    # Filtered dynamic spectrum
    dynspec = (dynspec3/spectrum).T

    return dynspec, dynspec1, spectrum, spect_raw, subband_frequencies, starttime, channel_timeresolution, channel_frequencies, channel_delays

#-------------------------------------------------------------------------------
# Input arguments
#-------------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument("-file", help="Input file", type=str, dest='file', nargs='?')
parser.add_argument("-dm", help="Dispersion measure", type=float, dest='dm', nargs='?', default=26.7641)
parser.add_argument("-npol", help="Number of polarizations", type=int, dest='npol', nargs='?', default=1)

args = parser.parse_args()

#-------------------------------------------------------------------------------
# Data extraction
#-------------------------------------------------------------------------------

if args.npol == 1 :
    nch   = 32
    t_int = 32
    dynspec, dynspec1, spectrum, spect_raw, fsb, t0, tch, fch, dch = data(outfilename0, DM=DM, nch=nch, t_int=t_int)
else :
    print("Two polarizations plot not available yet")

#-------------------------------------------------------------------------------
# Plotting
#-------------------------------------------------------------------------------

# Defining min, max color scale
cmin = min([np.median(dynspec)-np.std(dynspec), 0])
cmax = np.median(dynspec)+np.std(dynspec)

ch_freq = np.linspace(min(fch),max(fch),dynspec.shape[0])
for line in np.where(spectrum>2*np.median(spectrum))[0] :
    dynspec[line] = np.median(dynspec) * np.ones(dynspec.shape[1])
for line in np.where(spectrum<np.median(spectrum)/2)[0] :
    dynspec[line] = np.median(dynspec) * np.ones(dynspec.shape[1])

plt.figure(figsize=(10,7))
plt.imshow(dynspec,aspect='auto',vmin=cmin, vmax=cmax, origin='lower', interpolation='nearest', cmap='hot', extent=[0,dynspec.shape[1]*tch*t_int,min(fch),max(fch)])
plt.colorbar()

plt.savefig(file.replace('.h5', '.pdf'))
