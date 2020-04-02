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
import glob
from optparse import OptionParser


# Defining the workspace parameters
deg = math.pi / 180.
pi2 = math.pi / 2.

#------------------------------------------------------------------------
# Defining functions
#------------------------------------------------------------------------

def utc2jd(utctimestamp):
    from pycrtools import tools
    import datetime
    dtt=datetime.datetime.utcfromtimestamp(utctimestamp)
    timestring=datetime.datetime.strftime(dtt,'%Y-%m-%d %H:%M:%S')
    return tools.strdate2jd(timestring)

def findstem(arr):
    """Finding common substring to list of strings"""
    n = len(arr) 
    s = arr[0] 
    l = len(s) 
    res = "" 
    for i in range( l) : 
        for j in range( i + 1, l + 1) : 
            stem = s[i:j] 
            k = 1
            for k in range(1, n):  
                if stem not in arr[k]: 
                    break
            if (k + 1 == n and len(res) < len(stem)): 
                res = stem 
    return res 

def outfilenames(outfiledir, files, pol):
    """Defining output file names"""

    if len(files) > 1:
        fname = findstem(files).split('/')[-1] + 'tbb'
    else :
        fname = files[0].split('/')[-1][0:-3]
    outbf = outfiledir + fname + '_bf_pol'+str(pol)+'.h5'
    outdynspec = outfiledir + fname + '_fft_pol'+str(pol)
    return outbf, outdynspec

class beamformer:
    

    def __init__(self, infile, outfilename):
        
        # Input h5py data
        self.infile = infile #h5py.File(infile,'r')
        # Output filename
        self.outfilename = outfilename
        # Output h5py data. TODO: convert to pycrtools data
        self.outfile = h5py.File(self.outfilename,'w')     
        # Station
        self.station = ''
        self.station_name = ''
        self.station_number = 0

    def __station_selection(self):
        """Selecting station from inpyt file keywords."""
        # Select stations in file. We now only select the first station. 
        stations = [k for f in self.tbb_files for k in f.keys() if 'STATION' in k]

        # Checking that all stations are the same
        if not all(elem == stations[0] for elem in stations):
            print("Stations:\n", stations)
            sys.exit("ERROR: not all stations are the same")

        self.station = stations[0]
        self.station_name = self.station[-5:]
        self.station_number = int(self.station[-3:]) 

    def __header_dictionary(self):
        """Creating beamformed data header"""

        # TODO: use setHeader and updateHeader with the reqired keywords

    def __keyword_dictionary(self):
        """Creating output keyword dictionary"""

        s = self.station

        # Checking that all attributes are the same
        for k in self.tbb_files[0].attrs.keys() :
            if 'FILENAME' not in k :
                try:
                    if not all(f.attrs[k] == self.tbb_files[0].attrs[k] 
                           for f in self.tbb_files):
                        print("ERROR: Attribute ",k, " not matching every file")
                        sys.exit()
                except:
                    if not all(f.attrs[k][i] == att for f in self.tbb_files 
                           for i,att in enumerate(self.tbb_files[0].attrs[k])):
                        print("ERROR: Attribute ",k, " not matching every file")
                        sys.exit()

        # Keywords from input file
        # TODO: Check that all attributes are the same for each board
        for k in self.tbb_files[0].attrs.keys():
            try:
                self.outfile.attrs.create(k,data=self.tbb_files[0].attrs[k])
            except:
                if k in ["TARGETS","OBSERVATION_STATIONS_LIST"]:
                    self.outfile.attrs.create(k,
                        data=[vv for vv in self.tbb_files[0].attrs[k]])
                else:
                    print("failed to create key",k)

        # Additional keywords
        for k in ["DIPOLE_NAMES","SELECTED_DIPOLES"]:
            try :
                self.outfile.attrs.create(k,data=[str(d.replace('DIPOLE_', '')) 
                     for d in self.tbb_files[0][s].keys() if int(d[-3:])%2==pol])
            except :
                print("failed to create key",k)

        # Creating beamformed data groups
        self.outfile.create_group(s)
        self.outfile[s].create_group('WEIGHTS')
        self.outfile[s].create_group('BFDATA')
        for k in self.tbb_files[0][s].attrs.keys():
            try:
                if k == u'NOF_DIPOLES':
                    ndipoles = 0
                    for f in self.tbb_files :
                        ndipoles += f[s].attrs[k]
                    self.outfile[s].attrs.create(k,data=ndipoles)
                else :
                    self.outfile[s].attrs.create(k,
                        data=self.tbb_files[0][s].attrs[k])
            except:
                print("failed to create key",k)

    def read_calibration(
            self, inputfile, 
            caltabledir='/data/holography/Holog-20191212-1130/caltables/'
            ):
        """Reading calibration tables"""

        # TODO: check this
        filter_selection = inputfile.attrs[u'FILTER_SELECTION'].replace('A_','A-')
        caltablename = caltabledir + "CalTable-" + "" \
                       + "{:03d}".format(self.station_number) \
                       + '-' + filter_selection + '.dat'
        caltable=ct.readTable(caltablename)[1]
        
        return caltable        

    def __dipoles_subbands(self):
        """Function selecting dipoles and subbands from input data"""

        s = self.station
        # Selecting dipoles
        self.dipoles = [f[s].keys() for f in self.tbb_files]

        # Select subbands
        self.subbands = np.unique(np.concatenate([f[s][d].keys() 
                        for f in self.tbb_files for d in f[s].keys()]))
        self.subbands.sort()

        #return dipoles, subbands

    def __subband_offsets(self, sb):
        """Computing time offsets for each dipole of the subband"""

        print("Subband ", sb)
        s = self.station

        # Deleting previously exixting variables
        try :
            del available_dipoles, d0, num, timeresolution, starttimes 
            del datalengths, minstarttime, maxstarttime, diffstarttimes
            del offsets, datalength
        except:
            print("")

        # Select all dipoles that have this subband in their keys, 
        # for the even or odd polarisations
        available_dipoles = [d for f in self.tbb_files for d in f[s].keys() 
                            if sb in f[s][d] if int(d[-3:])%2==pol]
        available_dipoles.sort()

        # Calculating offsets
        # TODO: replace self.tbb_files[0]
        d0 = available_dipoles[0]
        #timeresoilution = self.tbb_files[0][s][d0][sb].attrs['TIME_RESOLUTION']
        num = np.where(np.char.find(self.dipoles, d0)!= -1)[0][0]
        timeresolution = self.tbb_files[num][s][d0][sb].attrs['TIME_RESOLUTION']
        starttimes = [(f[s][d][sb].attrs['TIME'], 
                     timeresolution*f[s][d][sb].attrs['SLICE_NUMBER']) 
                     for i,f in enumerate(self.tbb_files)
                     for d in self.dipoles[i] if d in available_dipoles]
        datalengths = [f[s][d][sb].shape[0] 
                      for i,f in enumerate(self.tbb_files)
                      for d in self.dipoles[i] if d in available_dipoles]

        minstarttime = min(starttimes)
        maxstarttime = max(starttimes)
        diffstarttimes = (maxstarttime[0] - minstarttime[0] 
                         + maxstarttime[1] - minstarttime[1]) / timeresolution

        offsets2  = [int(math.ceil(((st[0] - minstarttime[0]) + 
                    (st[1]-minstarttime[1])) / timeresolution)) 
                    for st in starttimes]
        flag_offsets = [num for num,o in enumerate(offsets2) 
                       if o>offset_max_allowed]
        available_dipoles = [d for num,d in enumerate(available_dipoles) 
                            if num not in flag_offsets]
        starttimes = [(f[s][d][sb].attrs['TIME'], 
                     timeresolution * f[s][d][sb].attrs['SLICE_NUMBER']) 
                     for i,f in enumerate(self.tbb_files) 
                     for d in self.dipoles[i] 
                     if d in available_dipoles]
        minstarttime = min(starttimes)
        maxstarttime = max(starttimes)

        offsets = [int(math.ceil((maxstarttime[0] - st[0] 
                  + maxstarttime[1] - st[1])/timeresolution)) 
                  for st in starttimes]
        datalength = min(datalengths) - max(offsets)

        print("subband",sb,"# Available dipoles",len(available_dipoles))

        return offsets, available_dipoles, datalength

    def __subband_beamforming(
            self, sb, offsets, available_dipoles, datalength, caltable
            ):
        """Beamforming one subband"""

        s = self.station
        sbnr = int(sb[-3:])

        print("Creating output datasets")
        print(sb, datalength)
        bfdata = np.zeros(shape = (datalength,), dtype = np.complex)
        sbdata_complex = np.zeros(shape = (datalength,), dtype=np.complex)
        d0 = available_dipoles[0]
        num = np.where(np.char.find(self.dipoles, d0)!= -1)[0][0]

        # Calculate time for the weights. 
        t = utc2jd(self.tbb_files[num][s][d0][sb].attrs['TIME'])
        weights = azel2beamweights(getRADEC_2_AZEL(ra,dec,t), 
                  self.station_name, 
                  self.tbb_files[num][s][d0][sb].attrs[u'CENTRAL_FREQUENCY']
                  ).toNumpy()
        ndipoles = len(available_dipoles) 

        # The beamforming part
        sbweight = np.zeros(shape = datalength)
        for d,o in zip(available_dipoles, offsets):
            board = np.where(np.char.find(self.dipoles, d) != -1)[0][0]
            f = self.tbb_files[board]
            if d in ['DIPOLE_145000000','DIPOLE_145000002']:
                continue
            print("Analysing dipole",d,"at offset",o,)
            dnr = int(d[-3:])

            # Read data from the offset, so that all are equal
            sbdata = f[s][d][sb][o:o+datalength]

            # Convert to complex numbers
            try:
                sbdata_complex.real = sbdata['real']
                sbdata_complex.imag = sbdata['imag']
            except:
                sbdata_complex = sbdata

            # In the beamforming we correct for the calibration delay and the beamforming weights
            # We're weighing by the number of available dipoles, as this can be different per dipole

            dipoledata = sbdata_complex * caltable[dnr,sbnr] * weights[dnr]

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
        """Adding subband metadata"""

        s = self.station
        d0 = available_dipoles[0]
        num = np.where(np.char.find(self.dipoles, d0)!= -1)[0][0]  
        
        # TODO: check that all files have same attributes
        for k in list(self.tbb_files[num][s][d0][sb].attrs):
            if k not in ['FLAG_OFFSETS','SLICE_NUMBER']:
                self.outfile[s]['BFDATA'][sb].attrs.create(k,
                    data=self.tbb_files[num][s][d0][sb].attrs[k])

        # Correct slice number for the offsets
        self.outfile[s]['BFDATA'][sb].attrs.create('SLICE_NUMBER', 
             data = self.tbb_files[num][s][d0][sb].attrs['SLICE_NUMBER'] 
                    + np.array(offsets)[np.array(available_dipoles)==d0][0])
        for k in [u'STATION_ID', u'NYQUIST_ZONE', u'SAMPLE_FREQUENCY']:
             self.outfile[s]['BFDATA'][sb].attrs.create(k,
                 data=self.tbb_files[num][s][d0].attrs[k])

        for k in [u'RCU_ID', u'RSP_ID']:
            self.outfile[s]['BFDATA'][sb].attrs.create(k,
                data=[f[s][d].attrs[k] for i,f in enumerate(self.tbb_files) 
                for d in available_dipoles if d in self.dipoles[i]])

        for k in [u'DIPOLE_IDS']:
            self.outfile[s]['BFDATA'][sb].attrs.create(k,data=[str(d) 
                for d in available_dipoles])
        

    def beamforming(self):
        """
        Beamforming function.
        It will generate an .h5 file with the tbb beamformed data.
        """

        # Opening files
        self.tbb_files = []
        for f in self.infile :
            self.tbb_files.append(h5py.File(f, 'r'))

        # Defining stations
        self.__station_selection() 
        s = self.station
        station_name = self.station_name
        station_number = self.station_number

        # Creating dictionary fo output data
        # todo: Fix __keyword_dictionary function to take several input files
        self.__keyword_dictionary()

        # Read calibration tables. 
        # These should also be in a directory in the image, but can be old.
        # TODO: check that it's the same for each file
        caltable = self.read_calibration(self.tbb_files[0])

        # Selecting dipoles and subbands
        self.__dipoles_subbands()

        # TODO loop over subbands, first a porton, then all
        subband = self.subbands[0]
        for sb in self.subbands[0:]:
            # Select all dipoles that have this subband in their keys, 
            # for the even or odd polarisations
            available_dipoles = [d for f in self.tbb_files for d in f[s].keys() 
                                if sb in f[s][d].keys() if int(d[-3:])%2==pol]
            available_dipoles.sort()
            if len(available_dipoles)==0:
                continue

            # Computing offsets
            offsets, available_dipoles, datalength = self.__subband_offsets(sb)
            
            # BEAMFORMING
            self.__subband_beamforming(sb, offsets, available_dipoles, 
                datalength, caltable)

            # Adding metadata to output file
            self.__subband_metadata(sb, offsets, available_dipoles) 

        # Closing files
        for f in self.tbb_files :
            f.close()

        self.outfile.close()

    def convert2beam(self, dynspecfile, DM=26.8, nch=16, t_int=6):
        """
        Creating dynamic spectrum from tbb beamformed data
        INPUT:
        dynspecfile : dynamic spectrum output file
        DM          : dispersion measure of the source
        nch         : number of channels
        t_int       : integration time
        """

        # TODO: FIX THIS FUNCTION

	# Opening beamformed data
        self.outfile = h5py.File(self.outfilename,'r')

        print(self.outfile.keys())
        st = list(self.outfile.keys())[0]

        subbands = sorted(list(self.outfile[st]['BFDATA'].keys()))
        nsb = len(subbands)
        datalen = self.outfile[st]['BFDATA'][subbands[0]].shape[0]
        datalen = datalen-datalen%(2*nch*t_int)


        complexdynspec = np.zeros(
                         dtype = self.outfile[st]['BFDATA'][subbands[0]].dtype, 
                         shape = (nsb,datalen))
        self.subband_frequencies = np.zeros(dtype = float,shape=(nsb))
        self.weights = np.zeros(
                       dtype = self.outfile[st][u'WEIGHTS'][subbands[0]].dtype, 
                       shape = (nsb,datalen))

        for num,sb in enumerate(subbands):
            complexdynspec[num] = self.outfile[st]['BFDATA'][sb][0:datalen]
            self.subband_frequencies[num] = (self
                .outfile[st]['BFDATA'][sb].attrs[u'CENTRAL_FREQUENCY'])
            #weights[num] = self.outfile[st][u'WEIGHTS'][sb][0:datalen]

        # Same resolution as BF data : nch=16
        data01 = complexdynspec[:,:].T
        s = data01.shape
        fftdata01 = np.fft.fft(data01.reshape(s[0]//nch, nch, 
                    s[1]).swapaxes(1,2), axis=2).reshape(s[0]//nch, nch*s[1])

        # Dynamic spectrum
        self.subband_bandwidth = (self
            .outfile[st]['BFDATA'][subbands[0]].attrs['BANDWIDTH']) #Hz
        self.subband_timeresolution = (self
            .outfile[st]['BFDATA'][subbands[0]].attrs['TIME_RESOLUTION']) #s

        self.channel_timeresolution = self.subband_timeresolution*nch
        channel_bandwidth = self.subband_bandwidth/nch
        self.channel_frequencies = np.zeros(dtype=float,shape=(nsb,nch))
        self.channel_delays = np.zeros(dtype=int,shape=(nsb,nch))

        for num,sbfreq in enumerate(self.subband_frequencies):
            self.channel_frequencies[num] = (np.arange(nch)-nch/2) \
                * channel_bandwidth + sbfreq
            self.channel_delays[num] = dd.freq_to_delay(DM, 
                self.channel_frequencies[num], self.channel_timeresolution)
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
        self.starttime = [self.outfile[st]['BFDATA'][sb].attrs['TIME']%100 
                         + self.outfile[st]['BFDATA'][sb].attrs['SLICE_NUMBER'] 
                         * self.subband_timeresolution for sb in sorted(subbands)]

        # Downsampling data
        f_int = nch//2 # Best visualized when nch/f_int=2
        dynspec2 = np.sum(dynspec1.reshape(dynspec1.shape[0]//t_int, 
                   t_int, dynspec1.shape[1])/t_int, axis=1)
        dynspec3 = np.sum(dynspec2.reshape(dynspec2.shape[0], 
                   dynspec2.shape[1]//f_int, f_int)/f_int, axis=2)

        self.spectrum = np.average(dynspec3,axis=0)

        # Writing beamformed data
        #self.dynspec = cr.hArray(float, name=dynspecfile, fill=(dynspec3/self.spectrum).T)
        #self.dynspec.write(dynspecfile)
        
        #self.dynspec = (dynspec3/self.spectrum).T
        #np.save(dynspecfile.replace('.h5', '.npy'), self.dynspec)

        beam = cr.hArray(fftdata01, ext='.beam')
	beam.write(dynspecfile)
        # TODO: Use here __header_keyword function
        #beam.writeheader()
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

# Input file
files=glob.glob(args[0])
files.sort()
files.insert(0, files.pop(-1))
print('Opening following files:')
for f in files : print(f)

# Checking if file exists or is a list
assert len(files) >= 1, "Input file(s) not found"

# Other inputs
(ra,dec) = (options.ra, options.dec)
pol = options.pol
assert(pol in [0,1])

offset_max_allowed = options.offset_max_allowed

# Defining output file names
outfilename, outfiledyn = outfilenames(options.outfiledir, files, pol)

#------------------------------------------------------------------------
# Beamforming
#------------------------------------------------------------------------

b = beamformer(infile = files, outfilename = outfilename)
b.beamforming()
b.convert2beam(outfiledyn)


print("File written: ", outfilename)
print("File written: ", outfiledyn)
