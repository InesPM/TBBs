from __future__ import print_function
from __future__ import division
import h5py
import pycrtools as cr
import numpy as np
import caltable as ct
from radec2azel import *
import dedispersion as dd
import math
import sys
from optparse import OptionParser

#------------------------------------------------------------------------
# Defining functions
#------------------------------------------------------------------------

def utc2jd(utctimestamp):
    from pycrtools import tools
    import datetime
    dtt=datetime.datetime.utcfromtimestamp(utctimestamp)
    timestring=datetime.datetime.strftime(dtt,'%Y-%m-%d %H:%M:%S')
    return tools.strdate2jd(timestring)

class beamformer:
    
    def __init__(self, infile, outfilename):
        
        # Input h5py data
        self.infile = infile
        # Output filename
        self.outfilename = outfilename
        # Output h5py data. TODO: convert to pycrtools data
        self.outfile = h5py.File(self.outfilename,'w')     
        # Station
        self.station = ''
        self.station_name = ''
        self.sttion_number = 0

    def __station_selection(self):
        """
        Selecting station from inpyt file keywords.
        """
        # Select stations in file. We now only select the first station. 
        # We could loop over stations, but at the moment one station is written per file
        stations = [k for k in self.infile.keys() if 'STATION' in k]
        self.station = stations[0]
        self.station_name = self.station[-5:]
        self.station_number = int(self.station[-3:]) 

    def __keyword_dictionary(self):
        """
        Creating output keyword dictionary
        """

        s = self.station

        # Keywords from input file
        for k in self.infile.attrs.keys():
            try:
                self.outfile.attrs.create(k,data=self.infile.attrs[k])
            except:
                if k in ["TARGETS","OBSERVATION_STATIONS_LIST"]:
                    self.outfile.attrs.create(k,data=[vv for vv in self.infile.attrs[k]])
                else:
                    print("failed to create key",k)

        # Additional keywords
        for k in ["DIPOLE_NAMES","SELECTED_DIPOLES"]:
            try :
                self.outfile.attrs.create(k,data=[str(d.replace('DIPOLE_', '')) 
                     for d in self.infile[station].keys() if int(d[-3:])%2==pol])
            except :
                print("failed to create key",k)

        # Creating beamformed data groups
        self.outfile.create_group(s)
        self.outfile[s].create_group('WEIGHTS')
        self.outfile[s].create_group('BFDATA')
        for k in self.infile[s].attrs.keys():
            try:
                self.outfile[s].attrs.create(k,data=self.infile[s].attrs[k])
            except:
                print("failed to create key",k)

    def read_calibration(self, caltabledir='/data/holography/Holog-20191212-1130/caltables/'):
        """
        Reading calibration tables
        """
        filter_selection = self.infile.attrs[u'FILTER_SELECTION'].replace('A_','A-')
        caltablename = caltabledir + "CalTable-" + "" + "{:03d}".format(self.station_number) \
                       + '-' + filter_selection + '.dat'
        caltable=ct.readTable(caltablename)[1]
        
        return caltable        

    def __dipoles_subbands(self):
        """
        Function selecting dipoles and subbands from input data
        """
        s = self.station
        # Selecting dipoles
        dipoles=self.infile[s].keys()
        
        # Select subbands
        subbands=np.unique(np.concatenate([self.infile[s][d].keys() for d in dipoles]))
        subbands.sort()

        return dipoles, subbands

    def __subband_offsets(self, sb, dipoles):
        """
        Computing time offsets for each dipole of the subband
        """
        s = self.station

        # Select all dipoles that have this subband in their keys, for the even or odd polarisations
        available_dipoles = [d for d in dipoles if sb in self.infile[s][d].keys() 
                            if int(d[-3:])%2==pol ]
        available_dipoles.sort()
        #if len(available_dipoles)==0:
        #    continue

        # Calculating offsets
        d = available_dipoles[0]
        timeresolution = self.infile[s][d][sb].attrs['TIME_RESOLUTION']
        starttimes = [(self.infile[s][d][sb].attrs['TIME'], 
                     timeresolution*self.infile[s][d][sb].attrs['SLICE_NUMBER']) 
                     for d in available_dipoles]
        datalengths = [self.infile[s][d][sb].shape[0] for d in available_dipoles]

        minstarttime = min(starttimes)
        maxstarttime = max(starttimes)
        diffstarttimes = (maxstarttime[0] - minstarttime[0] + maxstarttime[1] - minstarttime[1]) \
                         / timeresolution

        offsets2 = [int(math.ceil(((st[0] - minstarttime[0]) + 
                   (st[1]-minstarttime[1])) / timeresolution)) for st in starttimes]
        flag_offsets = [num for num,o in enumerate(offsets2) if o>offset_max_allowed]
        available_dipoles = [d for num,d in enumerate(available_dipoles) if num not in flag_offsets]
        starttimes = [(self.infile[s][d][sb].attrs['TIME'], 
                     timeresolution*self.infile[s][d][sb].attrs['SLICE_NUMBER']) 
                     for d in available_dipoles]
        minstarttime = min(starttimes)
        maxstarttime = max(starttimes)

        offsets = [int(round((maxstarttime[0] - st[0] + maxstarttime[1] - st[1])/timeresolution)) 
                  for st in starttimes]
        datalength = min(datalengths) - max(offsets)

        print("subband",sb,"# Available dipoles",len(available_dipoles))

        return offsets, available_dipoles, datalength

    def __subband_beamforming(self, sb, offsets, available_dipoles, datalength, caltable):
        """
        Beamforming one subband
        """

        s = self.station
        sbnr = int(sb[-3:])

        print("Creating output datasets")
        bfdata = np.zeros(shape = (datalength,), dtype = np.complex)
        sbdata_complex = np.zeros(shape = (datalength,), dtype=np.complex)
        d0 = available_dipoles[0]
        # Calculate time for the weights. 
        # We could vary the weight over the observation if this is longer than one second
        t = utc2jd(self.infile[s][d0][sb].attrs['TIME'])
        weights = azel2beamweights(getRADEC_2_AZEL(ra,dec,t), 
                  self.station_name, 
                  self.infile[s][d0][sb].attrs[u'CENTRAL_FREQUENCY']).toNumpy()
        ndipoles = len(available_dipoles) 

        # The beamforming part
        sbweight = np.zeros(shape=datalength)
        for d,o in zip(available_dipoles,offsets):
            if d in ['DIPOLE_145000000','DIPOLE_145000002']:
                continue
            print("Analysing dipole",d,"at offset",o,)
            dnr = int(d[-3:])

            # Read data from the offset, so that all are equal
            sbdata = self.infile[s][d][sb][o:o+datalength]

            # Convert to complex numbers
            try:
                sbdata_complex.real = sbdata['real']
                sbdata_complex.imag = sbdata['imag']
            except:
                sbdata_complex = sbdata

            # In the beamforming we correct for the calibration delay and the beamforming weights
            # We're weighing by the number of available dipoles, as this can be different per dipole

            dipoledata = sbdata_complex*caltable[dnr,sbnr]*weights[dnr]

            dipoleweight = dipoledata/dipoledata
            dipoleweight[np.isnan(dipoleweight)] = 0
            sbweight += dipoleweight
            nzeroes = np.sum(np.abs(dipoledata)==0)
            bfdata += dipoledata
            print("nzeroes",nzeroes/len(bfdata))

        sbweight[sbweight==0] = 1
        bfdata /= np.sqrt(sbweight)
        # Create subband dataset
        self.outfile[s]['BFDATA'].create_dataset(sb,data=bfdata)
        self.outfile[s]['WEIGHTS'].create_dataset(sb,data=sbweight)


    def __subband_metadata(self, sb, offsets, available_dipoles):
        """
        Adding subband metadata
        """

        s = self.station
        d = available_dipoles[0]
 
        for k in list(self.infile[s][d][sb].attrs):
            if k not in ['FLAG_OFFSETS','SLICE_NUMBER']:
                self.outfile[s]['BFDATA'][sb].attrs.create(k,data=self.infile[s][d][sb].attrs[k])

        # Correct slice number for the offsets
        self.outfile[s]['BFDATA'][sb].attrs.create('SLICE_NUMBER', 
             data=self.infile[s][d][sb].attrs['SLICE_NUMBER'] + 
             np.array(offsets)[np.array(available_dipoles)==d][0])
        for k in [
          u'STATION_ID',
          u'NYQUIST_ZONE',
          u'SAMPLE_FREQUENCY']:
             self.outfile[s]['BFDATA'][sb].attrs.create(k,data=self.infile[s][d].attrs[k])

        for k in [u'RCU_ID', u'RSP_ID']:
            self.outfile[s]['BFDATA'][sb].attrs.create(k,data=[self.infile[s][d].attrs[k] 
                 for d in available_dipoles])

        for k in [u'DIPOLE_IDS']:
            self.outfile[s]['BFDATA'][sb].attrs.create(k,data=[str(d) for d in available_dipoles])
        

    def beamforming(self):
        """
        Beamforming function
        """

        # Defining stations
        self.__station_selection() 
        s = self.station
        station_name = self.station_name
        station_number = self.station_number

        # Creating dictionary fo output data
        self.__keyword_dictionary()

        # Read calibration tables. These should also be in a directory in the image, but can be old.
        caltable = self.read_calibration()

        # Selecting dipoles and subbands
        dipoles, subbands = self.__dipoles_subbands()

        # TODO loop over subbands, first a porton, then all
        subband=subbands[0]
        for sb in subbands[0:2]:
            # Select all dipoles that have this subband in their keys, 
            # for the even or odd polarisations
            available_dipoles = [d for d in dipoles if sb in self.infile[s][d].keys() 
                                if int(d[-3:])%2==pol ]
            available_dipoles.sort()
            if len(available_dipoles)==0:
                continue

            # Computing offsets
            offsets, available_dipoles, datalength = self.__subband_offsets(sb, dipoles)
            
            # BEAMFORMING
            self.__subband_beamforming(sb, offsets, available_dipoles, datalength, caltable)

            # Adding metadata to output file
            self.__subband_metadata(sb, offsets, available_dipoles) 

        # Closing output file
        self.outfile.close()

    def convert2beam(self, dynspecfile, DM=26.8, nch=16, t_int=6):

	# Opening beamformed data
        self.outfile = h5py.File(self.outfilename,'r')

        print(self.outfile.keys())
        st = list(self.outfile.keys())[0]

        subbands = sorted(list(self.outfile[st]['BFDATA'].keys()))
        nsb = len(subbands)
        datalen = self.outfile[st]['BFDATA'][subbands[0]].shape[0]
        datalen = datalen-datalen%(2*nch*t_int)


        complexdynspec = np.zeros(dtype = self.outfile[st]['BFDATA'][subbands[0]].dtype, shape = (nsb,datalen))
        self.subband_frequencies = np.zeros(dtype = float,shape=(nsb))
        self.weights = np.zeros(dtype = self.outfile[st][u'WEIGHTS'][subbands[0]].dtype, shape = (nsb,datalen))

        for num,sb in enumerate(subbands):
            complexdynspec[num] = self.outfile[st]['BFDATA'][sb][0:datalen]
            self.subband_frequencies[num] = self.outfile[st]['BFDATA'][sb].attrs[u'CENTRAL_FREQUENCY']
            #weights[num] = self.outfile[st][u'WEIGHTS'][sb][0:datalen]

        # Same resolution as BF data : nch=16
        data01 = complexdynspec[:,:].T
        s = data01.shape
        fftdata01 = np.fft.fft(data01.reshape(s[0]//nch, nch, s[1]).swapaxes(1,2), axis=2).reshape(s[0]//nch, nch*s[1])

        # Dynamic spectrum
        self.subband_bandwidth = self.outfile[st]['BFDATA'][subbands[0]].attrs['BANDWIDTH'] #Hz
        self.subband_timeresolution = self.outfile[st]['BFDATA'][subbands[0]].attrs['TIME_RESOLUTION'] #seconds


        self.channel_timeresolution = self.subband_timeresolution*nch
        channel_bandwidth = self.subband_bandwidth/nch
        self.channel_frequencies = np.zeros(dtype=float,shape=(nsb,nch))
        self.channel_delays = np.zeros(dtype=int,shape=(nsb,nch))

        for num,sbfreq in enumerate(self.subband_frequencies):
            self.channel_frequencies[num] = (np.arange(nch)-nch/2) * channel_bandwidth + sbfreq
            self.channel_delays[num] = dd.freq_to_delay(DM, self.channel_frequencies[num], self.channel_timeresolution)
        self.channel_frequencies = self.channel_frequencies.reshape((nsb*nch))
        self.channel_delays = self.channel_delays.reshape((nsb*nch))

        dynspec1 = np.abs(fftdata01)
        if False:
            for ch,delay in enumerate(self.channel_delays):
                if int(delay) >= 1:
                    dynspec1[0:-int(delay),ch] = dynspec1[int(delay):,ch]
                    dynspec1[-int(delay):,ch] = 0 #dynspec[0:int(delay),ch]

        # Spectrum
        specstd = np.std(dynspec1,axis=0)
        spect_raw = np.sum(dynspec1,axis=0)
        noisyness = (spect_raw/np.median(spect_raw))/(specstd/np.median(specstd))
        flagfreq = (noisyness<0.95)|(noisyness>1.05)

        # Start time
        self.starttime = [self.outfile[st]['BFDATA'][sb].attrs['TIME']%100 + self.outfile[st]['BFDATA'][sb].attrs['SLICE_NUMBER'] * self.subband_timeresolution for sb in sorted(subbands)]

        # Downsampling data
        f_int = nch//2 # Best visualized when nch/f_int=2
        dynspec2 = np.sum(dynspec1.reshape(dynspec1.shape[0]//t_int, t_int, dynspec1.shape[1])/t_int, axis=1)
        dynspec3 = np.sum(dynspec2.reshape(dynspec2.shape[0], dynspec2.shape[1]//f_int, f_int)/f_int, axis=2)

        self.spectrum = np.average(dynspec3,axis=0)

        # Writing beamformed data
        #self.dynspec = cr.hArray(float, name=dynspecfile, fill=(dynspec3/self.spectrum).T)
        #self.dynspec.write(dynspecfile)
        
        #self.dynspec = (dynspec3/self.spectrum).T
        #np.save(dynspecfile.replace('.h5', '.npy'), self.dynspec)

        beam = cr.hArray(fftdata01, ext='.beam')
	beam.write(dynspecfile)
        # Closing beamformed data
        self.outfile.close()

#------------------------------------------------------------------------
# Command line options
#------------------------------------------------------------------------
parser = OptionParser()

parser.add_option("-r", "--right_ascension", dest="ra", type="float", default=0.92934186635, help="Right ascension of the source in degrees. Default RA: B0329")
parser.add_option("-d", "--declination", dest="dec", type="float", default=0.952579228492, help="Declination of the source. Default Dec: B0329")
parser.add_option("-p", "--polarisation", dest="pol", type="int", default=0, help="Polarisation to beamform.")
parser.add_option("-o", "--offset_max_allowed", dest="offset_max_allowed", type="int", default=400, help="Maximum offset between subbands in time bins.")
parser.add_option("-f", "--outfiledir", type="str", default="/data/projects/COM_ALERT/pipeline/analysis/marazuela/data/", help="Directory where the output files will be generated.")

(options, args) = parser.parse_args()

ffull=args[0]

(ra,dec) = (options.ra, options.dec)
pol = options.pol
assert(pol in [0,1])

offset_max_allowed = options.offset_max_allowed
outfilename = options.outfiledir + ffull.split('/')[-1][0:-3]+'bf_pol'+str(pol)+'.h5'
outfiledyn = options.outfiledir + ffull.split('/')[-1][0:-3]+'_fft_pol'+str(pol)

#------------------------------------------------------------------------
# Opening input file
#------------------------------------------------------------------------

infile = h5py.File(ffull,'r')

#------------------------------------------------------------------------
# Beamforming
#------------------------------------------------------------------------

b = beamformer(infile, outfilename)
b.beamforming()
b.convert2beam(outfiledyn)


print("File written: ", outfilename)
print("File written: ", outfiledyn)
