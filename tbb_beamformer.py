from __future__ import print_function
from __future__ import division
import h5py
import pycrtools as cr
from pycrtools import Beam_Tools as bt
import numpy as np
import caltable as ct
from radec2azel import *
import dedispersion as dd
from math import *
import sys
import shutil
import os
import glob
from optparse import OptionParser
from datetime import datetime

#------------------------------------------------------------------------
# Defining functions
#------------------------------------------------------------------------

def utc2jd(utctimestamp):
    from pycrtools import tools
    import datetime
    dtt=datetime.datetime.utcfromtimestamp(utctimestamp)
    timestring=datetime.datetime.strftime(dtt,'%Y-%m-%d %H:%M:%S')
    return tools.strdate2jd(timestring)

class BeamFormer:
    

    def __init__(self, infile, bffilename, 
            pol=0, offset_max_allowed=400, overwrite=True):
        
        # Input h5py data
        self.infile = infile #h5py.File(infile,'r')

        # Output filename
        self.bffilename = bffilename # TODO: Replace by bffile

        # Other inputs
        self.overwrite = overwrite
        self.pol = pol
        self.offset_max_allowed = offset_max_allowed

    #-------------------------------------------------------------------
    # Internal functions called by the beamformer

    def __station_selection(self):
        """Selecting station from input file keywords."""

        # Select stations in file. We now only select the first station. 
        stations = [k for f in self.tbb_files for k in f.keys() 
                if 'STATION' in k]

        # Checking that all stations are the same
        if not all(elem == stations[0] for elem in stations):
            print("Stations:\n", stations)
            sys.exit("ERROR: not all stations are the same")

        self.station = stations[0]
        self.station_name = self.station[-5:]
        self.station_number = int(self.station[-3:]) 

    def __data_groups(self):
        """Creating output data groups"""

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
        for k in self.tbb_files[0].attrs.keys():
            if k not in ["TARGETS", "FILENAME", "OBSERVATION_STATIONS_LIST"]:
              try:
                self.bffile.attrs.create(k,data=self.tbb_files[0].attrs[k])
              except:
                print("failed to create key",k)

        s = self.station

        # Creating beamformed data groups
        self.bffile.create_group(s)
        self.bffile[s].create_group('WEIGHTS')
        self.bffile[s].create_group('BFDATA')

    def read_calibration(self,  
            caltabledir='/data/holography/Holog-20191212-1130/caltables/'):
        """Reading calibration tables"""

        # TODO: check this
        
        filter_selection = (self.tbb_files[0].attrs[u'FILTER_SELECTION']
                .replace('A_','A-'))
        caltablename = (caltabledir + "CalTable-" + "" 
                + "{:03d}".format(self.station_number) 
                + '-' + filter_selection + '.dat')
        caltable=ct.readTable(caltablename)[1]
        
        return caltable        

    def __dipoles_subbands(self):
        """Function selecting dipoles and subbands from input data"""

        s = self.station
        # Selecting dipoles
        self.dipoles = [f[s].keys() for f in self.tbb_files]
        self.ndipoles = len([d for dip in self.dipoles for d in dip])

        # Select subbands
        self.subbands = np.unique(np.concatenate([f[s][d].keys() 
                        for f in self.tbb_files for d in f[s].keys()]))
        self.subbands.sort()
        self.nsubbands = len(self.subbands)

        #return dipoles, subbands

    def __subband_offsets(self, sb):
        """Computing time offsets for each dipole of the subband"""

        s = self.station

        # Select all dipoles with the given subband for the given polarisation
        available_dipoles = [d for f in self.tbb_files for d in f[s].keys() 
                if sb in f[s][d] if int(d[-3:])%2==self.pol]
        available_dipoles.sort()

        # Calculating offsets
        d0 = available_dipoles[0]
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

        offsets2  = [int(ceil(((st[0] - minstarttime[0]) + 
                (st[1]-minstarttime[1])) / timeresolution)) 
                for st in starttimes]
        flag_offsets = [num for num,o in enumerate(offsets2) 
                if o > self.offset_max_allowed]
        available_dipoles = [d for num,d in enumerate(available_dipoles) 
                if num not in flag_offsets]
        starttimes = [(f[s][d][sb].attrs['TIME'], 
                timeresolution * f[s][d][sb].attrs['SLICE_NUMBER']) 
                for i,f in enumerate(self.tbb_files) 
                for d in self.dipoles[i] 
                if d in available_dipoles]
        minstarttime = min(starttimes)
        maxstarttime = max(starttimes)

        offsets = [int(ceil((maxstarttime[0] - st[0] 
                + maxstarttime[1] - st[1])/timeresolution)) 
                for st in starttimes]
        datalength = min(datalengths) - max(offsets)

        return offsets, available_dipoles, datalength

    def __subband_beamforming(
            self, sb, offsets, available_dipoles, datalength, caltable
            ):
        """
        Beamforming one subband

        In the beamforming we correct for the calibration delay 
        and the beamforming weights.
        We're weighing by the number of available dipoles, 
        as this can be different per dipole

        Input:
        sb: subband id.
        offsets: given by __subband_offsets
        available dipoles: dipoles with data for the given polarization
        datalength: number of time bins
        caltable: calibration table given by read_calibration

        Output:
        Creation of subband dataset

        """

        s = self.station
        sbnr = int(sb[-3:])

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
            dnr = int(d[-3:])

            # Read data from the offset, so that all are equal
            sbdata = f[s][d][sb][o:o+datalength]

            # Convert to complex numbers
            try:
                sbdata_complex.real = sbdata['real']
                sbdata_complex.imag = sbdata['imag']
            except:
                sbdata_complex = sbdata

            dipoledata = sbdata_complex * caltable[dnr,sbnr] * weights[dnr]
            dipoleweight = dipoledata/dipoledata
            dipoleweight[np.isnan(dipoleweight)] = 0
            nzeroes = np.sum(np.abs(dipoledata)==0)

            sbweight += dipoleweight
            bfdata += dipoledata

        sbweight[sbweight==0] = 1
        bfdata /= np.sqrt(sbweight)

        # Create subband dataset
        self.bffile[s]['BFDATA'].create_dataset(sb,data=bfdata)
        self.bffile[s]['WEIGHTS'].create_dataset(sb,data=sbweight)


    def __subband_metadata(self, sb, offsets, available_dipoles):
        """Adding subband metadata"""

        s = self.station
        d0 = available_dipoles[0]
        num = np.where(np.char.find(self.dipoles, d0)!= -1)[0][0]  
        
        # TODO: check that all files have same attributes
        for k in list(self.tbb_files[num][s][d0][sb].attrs):
            if k not in ['FLAG_OFFSETS','SLICE_NUMBER']:
                self.bffile[s]['BFDATA'][sb].attrs.create(k,
                    data=self.tbb_files[num][s][d0][sb].attrs[k])

        # Correct slice number for the offsets
        self.bffile[s]['BFDATA'][sb].attrs.create('SLICE_NUMBER', 
                 data = self.tbb_files[num][s][d0][sb].attrs['SLICE_NUMBER'] 
                 + np.array(offsets)[np.array(available_dipoles)==d0][0])
        for k in [u'STATION_ID', u'NYQUIST_ZONE', u'SAMPLE_FREQUENCY']:
             self.bffile[s]['BFDATA'][sb].attrs.create(k,
                 data=self.tbb_files[num][s][d0].attrs[k])

        for k in [u'RCU_ID', u'RSP_ID']:
            self.bffile[s]['BFDATA'][sb].attrs.create(k,
                data=[f[s][d].attrs[k] for i,f in enumerate(self.tbb_files) 
                for d in available_dipoles if d in self.dipoles[i]])

        for k in [u'DIPOLE_IDS']:
            self.bffile[s]['BFDATA'][sb].attrs.create(k,data=[str(d) 
                for d in available_dipoles])    

    def __keyword_dictionary(self):
        """Creating output keyword dictionary"""

        f = self.tbb_files[0]
        d0 = f[self.station].keys()[0]
        sb0 = f[self.station][d0].keys()[0]
        fds = f[self.station][d0][sb0]
        # Creating keyword dictionary

        self.keyword_dict = {

            # Keywords replaced from input file
            'FILENAME': self.bffilename.split('/')[-1],

            # Keywords from data attributes
            'TIME_RESOLUTION': fds.attrs['TIME_RESOLUTION'],
            'BANDWIDTH_UNIT': fds.attrs['BANDWIDTH_UNIT'],
            'SAMPLES_PER_FRAME': fds.attrs['SAMPLES_PER_FRAME'], 
            'TIME_RESOLUTION_UNIT': fds.attrs['TIME_RESOLUTION_UNIT'],
            'DATA_LENGTH': fds.attrs['DATA_LENGTH'], 
            'TIME': fds.attrs['TIME'],
            'BANDWIDTH': fds.attrs['BANDWIDTH'],

            # TODO: Add keywords from original attribute list
            'TARGETS': [vv for vv in self.tbb_files[0].attrs['TARGETS']],
            'OBSERVATION_STATIONS_LIST': 
                    [vv for vv 
                    in self.tbb_files[0].attrs['OBSERVATION_STATIONS_LIST']],

            # Derived keywords
            'STATION': self.station,
            'STATION_NAME': self.station_name,
            'STATION_NUMBER': self.station_number,
            'NOF_SELECTED_DATASETS': self.ndipoles,
            'NOF_DIPOLE_DATASETS': self.ndipoles,
            'DIPOLE_NAMES': [str(d) for dip in self.dipoles for d in dip],
            'SELECTED_DIPOLES': 
                    [str(d.replace('DIPOLE_','')) for dip in self.dipoles 
                    for d in dip],
            'SELECTED_DIPOLES_INDEX': range(self.ndipoles),
            'CHANNEL_ID': 
                    [int(d.replace('DIPOLE_','')) for dip in self.dipoles
                    for d in dip],
            'SUBBANDS': [str(sb) for sb in self.subbands],
            'NOF_SUBBANDS': self.nsubbands,
            'SAMPLE_FREQUENCY': 200000000.0,
            'SAMPLE_FREQUENCY_UNIT': f.attrs["CLOCK_FREQUENCY_UNIT"], 
            'SAMPLE_FREQUENCY_VALUE': 200.0,
            'BANDWIDTH': 
                    self.bffile[self.station]['BFDATA'][self.subbands[0]]
                    .attrs['BANDWIDTH'], # Check this
            #'CLOCK_FREQUENCY_UNIT': ,
            'SAMPLE_INTERVAL': 5e-9,
            'SAMPLE_INTERVAL_UNIT': 's',
            'DIPOLE_CALIBRATION_DELAY_UNIT': 's',
            'OBSERVER': "I. Pastor-Marazuela"
        }
       
        # Writing keyword dictionary to output file

        # Derived keywords
        for k in self.keyword_dict.keys():
            try:
              if isinstance(self.keyword_dict[k], unicode):
                self.bffile.attrs.create(k, data=str(self.keyword_dict[k]))
              elif isinstance(self.keyword_dict[k], list):
                self.bffile.attrs.create(k, 
                    data=[vv for vv in self.keyword_dict[k]])
              else:
                self.bffile.attrs.create(k, data=self.keyword_dict[k])
            except:
              print(self.keyword_dict[k], type(self.keyword_dict[k]))
              print("failed to create key",k)

    #-------------------------------------------------------------------
    # Beamforming data

    def beamforming(self, test=False):
        """
        Beamforming function.
        It will generate an .h5 file with the tbb beamformed data.
        """

        print("--------------------------")
        print("       BEAMFORMING        \n")

        # Opening output h5py data 
        if self.overwrite:
            if os.path.exists(self.bffilename):
                os.remove(self.bffilename)
        self.bffile = h5py.File(self.bffilename,'w')

        # Opening input data
        self.tbb_files = []
        for f in self.infile :
            self.tbb_files.append(h5py.File(f, 'r'))

        # Defining stations
        self.__station_selection() 
        s = self.station
        station_name = self.station_name
        station_number = self.station_number

        # Creating output data groups
        self.__data_groups()

        # Read calibration tables. 
        caltable = self.read_calibration()

        # Selecting dipoles and subbands
        self.__dipoles_subbands()

        # Defining subbands to loop over. Less subbands for a test
        if test :
            sbs = self.subbands[0:2]
        else :
            sbs = self.subbands[0:]

        # Loop over subbands
        for sb in sbs:

            print ("Subband {}/{}".format(sb, sbs[-1]))
            sys.stdout.write("\033[F") # Cursor up one line

            # Select all dipoles with the given subband and given polarisation
            available_dipoles = [d for f in self.tbb_files for d in f[s].keys()
                    if sb in f[s][d].keys() if int(d[-3:])%2==self.pol]
            available_dipoles.sort()
            if len(available_dipoles)==0:
                continue

            # Computing offsets
            offsets, available_dipoles, datalength = self.__subband_offsets(sb)
            # Beamforming
            self.__subband_beamforming(sb, offsets, available_dipoles, 
                    datalength, caltable)
            # Adding metadata to output file
            self.__subband_metadata(sb, offsets, available_dipoles)

        # Creating keyword dictionary 
        self.__keyword_dictionary() 

        # Closing files
        for f in self.tbb_files :
            f.close()
        self.bffile.close()

    #-------------------------------------------------------------------
    # FFT related functions

    def dedispersed_time(self, t0, f1, f2=0.2):
        """
        Getting dedispersed time for given time and frequency subband
        """
        f1 = f1*1e-9 # GHz
        t = 4.15e-3 * self.dm * (1/f1**2 - 1/f2**2)

        return t0-t

    def __get_time_hr(self):
        """
        Getting the 'TIME_HR' header keyword from the input data.
        Dedispersing to the reference frequency 200 MHz
        """

        s = self.station
        self.tbb_files = []
        for f in self.infile :
            self.tbb_files.append(h5py.File(f, 'r'))

        self.time_hr = []
        self.time_key = []
        self.sample_number = []

        for f in self.tbb_files:
            for d in f[s].keys():
                sbs = f[s][d].keys()
                tr = f[s][d][sbs[-1]].attrs['TIME_RESOLUTION']
                t0 = (f[s][d][sbs[-1]].attrs['TIME'] 
                        + tr * f[s][d][sbs[-1]].attrs['SLICE_NUMBER'])
                f1 = f[s][d][sbs[-1]].attrs[u'CENTRAL_FREQUENCY']

                t = self.dedispersed_time(t0, f1)
                self.time_key.append(int(t))

                tstr = datetime.fromtimestamp(t).strftime("%Y-%m-%d %H:%M:%S")
                self.time_hr.append(tstr)

                samplenum = int((t - floor(t))/5e-9)
                self.sample_number.append(samplenum)

    def __bf_dictionary(self):
        """
        Defining BeamFormer dictionary, 
        with keywords from pycrtools/beamformer.py
        required for imaging
        """

        stride = 1
        delta_nu = 1 
        maxchunklen = 2*20
        filesize = 204800000
        sample_interval = 5e-9
        blocklen = len(self.channel_frequencies)
        chunklen = min(filesize/stride, maxchunklen)
        max_nblocks = int(floor(filesize / stride / blocklen))
        #nblocks = int(min(max(round(chunklen / blocklen), 1), max_nblocks))
        speclen = len(self.channel_frequencies) #blocklen/2 + 1
        start_time = 0
        end_time = self.fftdata.shape[0]*5.12e-6 * self.nch
        block_duration = 5.12e-6 * self.nch
        nblocks = self.fftdata.shape[0]

        self.bf_dict = {
            'stride': stride,
            'delta_nu': delta_nu,
            'filesize': filesize,
            'sample_interval': sample_interval,
            'maxchunklen': filesize,
            'chunklen': chunklen,
            'max_nblocks': max_nblocks,

            'antenna_set': 'HBA1',
            'blocklen': blocklen,
            'nantennas_total': 48,
            'nblocks': nblocks,
            'speclen': speclen,

            'start_time': start_time,
            'end_time': end_time,
            'block_duration': block_duration
        }

    def write_dynspec(self) :
        """
        Saving Beamformed FFT data into a .beam file
        """

        # Writing header
        #f = h5py.File(self.infile[0], 'r')
        f = self.bffile
        ndipoles = f.attrs["NOF_DIPOLE_DATASETS"]

        self.__get_time_hr()
        self.__bf_dictionary()

        # Writing data        
        hfftdata = cr.hArray(self.fftdata, name="Beamed FFT")

        hfftdata.setHeader(

            # BeamFormer keywords
            BeamFormer = self.bf_dict,

            # Existing keywords
            OBSERVATION_START_UTC = str(f.attrs["OBSERVATION_START_UTC"]),
            OBSERVATION_ID = str(f.attrs["OBSERVATION_ID"]),
            CLOCK_FREQUENCY_UNIT = str(f.attrs["CLOCK_FREQUENCY_UNIT"]),
            NOTES = str(f.attrs["NOTES"]),
            OBSERVATION_FREQUENCY_CENTER =
                    f.attrs["OBSERVATION_FREQUENCY_CENTER"],
            PROJECT_PI = str(f.attrs["PROJECT_PI"]),
            OBSERVATION_END_UTC = str(f.attrs["OBSERVATION_END_UTC"]),
            PROJECT_CO_I = str(f.attrs["PROJECT_CO_I"]),
            TELESCOPE = str(f.attrs["TELESCOPE"]),
            ANTENNA_SET = str(f.attrs["ANTENNA_SET"]),
            OBSERVATION_START_MJD = f.attrs["OBSERVATION_START_MJD"],
            PROJECT_CONTACT = str(f.attrs["PROJECT_CONTACT"]),
            FILTER_SELECTION = str(f.attrs["FILTER_SELECTION"]),
            FILETYPE = str(f.attrs["FILETYPE"]),
            OBSERVATION_FREQUENCY_MAX = f.attrs["OBSERVATION_FREQUENCY_MAX"],
            CLOCK_FREQUENCY = f.attrs["CLOCK_FREQUENCY"],
            OBSERVATION_END_MJD = f.attrs["OBSERVATION_END_MJD"],
            OBSERVATION_NOF_STATIONS = f.attrs["OBSERVATION_NOF_STATIONS"],
            OBSERVATION_FREQUENCY_UNIT =
                    str(f.attrs["OBSERVATION_FREQUENCY_UNIT"]),
            SYSTEM_VERSION = str(f.attrs["SYSTEM_VERSION"]),
            OBSERVATION_FREQUENCY_MIN = f.attrs["OBSERVATION_FREQUENCY_MIN"],
            PROJECT_ID = str(f.attrs["PROJECT_ID"]),
            PROJECT_TITLE = str(f.attrs["PROJECT_TITLE"]),
            FILEDATE = str(f.attrs["FILEDATE"]),
            FILENAME = self.dynspecfile.split('/')[-1],
            TARGET = f.attrs["TARGETS"],

            # Keywords derived in the beamforming part
            STATION_NAME = [str(f.attrs["STATION_NAME"])] * ndipoles,
            DIPOLE_NAMES = f.attrs["DIPOLE_NAMES"],
            NOF_SELECTED_DATASETS = f.attrs["NOF_SELECTED_DATASETS"],
            NOF_DIPOLE_DATASETS = f.attrs["NOF_DIPOLE_DATASETS"],
            SELECTED_DIPOLES = f.attrs["SELECTED_DIPOLES"],
            SELECTED_DIPOLES_INDEX = range(ndipoles),
            CHANNEL_ID = f.attrs["CHANNEL_ID"],
            SAMPLE_FREQUENCY = [f.attrs["SAMPLE_FREQUENCY"]] * ndipoles,
            SAMPLE_FREQUENCY_UNIT = 
                    [str(f.attrs["SAMPLE_FREQUENCY_UNIT"])] * ndipoles,
            SAMPLE_FREQUENCY_VALUE = 
                    [f.attrs["SAMPLE_FREQUENCY_VALUE"]] * ndipoles,
            FREQUENCY_INTERVAL = [f.attrs["BANDWIDTH"]/self.nch] * ndipoles,
            SAMPLE_INTERVAL = [f.attrs["SAMPLE_INTERVAL"]] * ndipoles,
            SAMPLE_INTERVAL_UNIT = 
                    [str(f.attrs["SAMPLE_INTERVAL_UNIT"])] * ndipoles,
            DIPOLE_CALIBRATION_DELAY_UNIT = 
                    [str(f.attrs["DIPOLE_CALIBRATION_DELAY_UNIT"])] * ndipoles,
            OBSERVER = str(f.attrs["OBSERVER"]),

            # Derived keywords that we need
            FREQUENCY_DATA = 
                    cr.hArray(self.channel_frequencies, name="Frequency"),
            BEAM_FREQUENCIES = 
                    cr.hArray(self.channel_frequencies, name="Frequency"),
            ALIGNMENT_REFERENCE_ANTENNA = 0, # TODO: check how this is defined
            PIPELINE_NAME = "UNDEFINED",
            BLOCK = 0,
            BLOCKSIZE = "", # 1024
            clock_offset = 
                    [md.getClockCorrection(self.station_name,
                    antennaset='HBA_DUAL', time=t)
                    for t in self.time_key],
            MAXIMUM_READ_LENGTH = 204800000, # Check
            FFTSIZE = self.fftdata.shape[-1],
            SAMPLE_NUMBER = [177325056] * ndipoles, # Check 177325056
            NYQUIST_ZONE = [2] * ndipoles,
            FREQUENCY_RANGE = "", # [(100000000.0, 200000000.0)] * 96
            PIPELINE_VERSION = "UNDEFINED",

            # Derived keywords that are less relevant
            TIMESERIES_DATA = "", # Check
            #"EMPTY_FFT_DATA = cr.hArray(Type=complex, name="fft(E-Field)")
            #"CABLE_ATTENUATION =
            TIME_HR = self.time_hr,
            TIME = self.time_key,
            ANTENNA_POSITION = "",
                    #  [3826575.52551, 460961.8472, 5064899.465]*96
            NOF_STATION_GROUPS = 1, # Check
            ITRF_ANTENNA_POSITION = "",
                    # Same as ANTENNA_POSITION" but cr.hArray
            CABLE_LENGTH = cr.hArray([115] * ndipoles), # Check
            OBSERVATION_END_TAI = "UNDEFINED",
            OBSERVATION_STATION_LIST = ["UNDEFINED"],
            STATION_GAIN_CALIBRATION = "", # No idea
            LCR_ANTENNA_POSITION = "", # No idea
            CABLE_DELAY = "", # Error
            SCR_ANTENNA_POSITION = "", # No idea
            OBSERVATION_START_TAI = "UNDEFINED",
            DATA_LENGTH = [204800000] * ndipoles, # Check
            CABLE_DELAY_UNIT = "", # Error

            # Possibly generated keywords
            TIME_DATA = "", # No idea
            FFT_DATA = cr.hArray(self.fftdata, name="fft(E-Field)")
            #"EMPTY_TIMESERIES_DATA = "",

        )


        # Writing to file
        hfftdata.write(self.dynspecfile, writeheader=True, 
                clearfile=True, ext='beam')

        # Writing header
        hfftdata.writeheader(self.dynspecfile, ext='beam')

        #t = cr.open(self.dynspecfile)
        #print(t.keys())

    #-------------------------------------------------------------------
    # FFTing data and saving to .beam file

    def convert2beam(self, dynspecfile, DM=26.8, nch=16, t_int=6):
        """
        Creating dynamic spectrum from tbb beamformed data
        INPUT:
        dynspecfile : dynamic spectrum output file
        DM          : dispersion measure of the source
        nch         : number of channels
        t_int       : integration time
        """

        print("--------------------------")
        print("     CONVERT TO .BEAM     \n")

	# Opening beamformed data
        print("Opening ", self.bffilename)
        self.dynspecfile = dynspecfile
        self.bffile = h5py.File(self.bffilename,'r')

        if self.overwrite:
            if os.path.exists(self.dynspecfile):
                shutil.rmtree(self.dynspecfile)

        # Setting parameters
        self.dm = DM
        self.nch = nch
        self.t_int = t_int
        self.station = self.bffile.attrs["STATION"]
        self.station_name = self.bffile.attrs["STATION_NAME"]
        self.station_number = self.bffile.attrs["STATION_NUMBER"]
        self.dipoles = self.bffile.attrs["DIPOLE_NAMES"]
        self.subbands = self.bffile.attrs["SUBBANDS"]

        st = list(self.bffile.keys())[0]

        #self.subbands = sorted(list(self.bffile[st]['BFDATA'].keys()))
        nsb = len(self.subbands)
        sb0 = self.subbands[0]

        datalen = self.bffile[st]['BFDATA'][sb0].shape[0]
        datalen = datalen-datalen%(2*nch*t_int)

        complexdynspec = np.zeros(
                dtype = (self.bffile[st]['BFDATA'][sb0].dtype), 
                shape = (nsb,datalen))
        self.subband_frequencies = np.zeros(dtype = float,shape=(nsb))
        self.weights = np.zeros(
                dtype = self.bffile[st][u'WEIGHTS'][sb0].dtype, 
                shape = (nsb,datalen))

        for num,sb in enumerate(self.subbands):
            complexdynspec[num] = self.bffile[st]['BFDATA'][sb][0:datalen]
            self.subband_frequencies[num] = (self
                    .bffile[st]['BFDATA'][sb].attrs[u'CENTRAL_FREQUENCY'])
            #weights[num] = self.outfile[st][u'WEIGHTS'][sb][0:datalen]

        # Same resolution as BF data : nch=16
        data01 = complexdynspec[:,:].T
        sh = data01.shape
        self.fftdata = np.fft.fft(data01.reshape(sh[0]//nch, nch, 
                sh[1]).swapaxes(1,2), axis=2).reshape(sh[0]//nch, nch*sh[1])
        self.fftdata = self.fftdata.reshape(self.fftdata.shape[0], 1, 
                self.fftdata.shape[1])

        # Dynamic spectrum
        # There is a keyword for bandwidth in bffile. Also for time resolution
        self.subband_bandwidth = (self
                .bffile[st]['BFDATA'][sb0].attrs['BANDWIDTH']) #Hz
        self.subband_timeresolution = (self
                .bffile[st]['BFDATA'][sb0].attrs['TIME_RESOLUTION']) #s

        self.channel_timeresolution = self.subband_timeresolution * nch
        self.channel_bandwidth = self.subband_bandwidth / nch
        self.channel_frequencies = np.zeros(dtype=float,shape=(nsb,nch))
        self.channel_delays = np.zeros(dtype=int,shape=(nsb,nch))

        for num,sbfreq in enumerate(self.subband_frequencies):
            self.channel_frequencies[num] = ((np.arange(nch)-nch/2)
                    * self.channel_bandwidth + sbfreq)
            self.channel_delays[num] = dd.freq_to_delay(DM, 
                    self.channel_frequencies[num], self.channel_timeresolution)
        self.channel_frequencies = self.channel_frequencies.reshape((nsb*nch))
        self.channel_delays = self.channel_delays.reshape((nsb*nch))

        # Writing data to file
        self.write_dynspec()

        # Closing beamformed data
        self.bffile.close()

#------------------------------------------------------------------------
# End of the script
