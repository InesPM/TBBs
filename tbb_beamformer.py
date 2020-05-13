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
    """
    Description
    -----------
    Beamforming TBB data for one station.
    
    Main functions
    --------------

    beamforming:        
        Beamforming TBB data, output in .h5 format

    convert2beam:       
        FFTing beamformed data, output in .beam format that can be 
        beamformed across stations with the pycrtools addBeams function.
        To use once the beamformed .h5 data exists.
    """

    def __init__(self, infile, bffilename, ra, dec,
            pol=0, substation='HBA', offset_max_allowed=400, 
            lofar_centered=True, reftime=None, dm=0.0, overwrite=True, 
            test=False):
        """
        Parameters
        ----------
        infile : str             
            List of .h5 files with TBB data of one station.
        bffilename : str
            Name of the output beamformed .h5 file.
        ra : float
            Right ascension of the source.
        dec : float
            Declination of the source.
        pol : int (0|1)
            Polarization to beamform.
        substation : str (HBA0|HBA1|HBA)
            Substation to beamform. 
            (HBA0|HBA1) for core stations. HBA for remote stations.
        offset_max_allowed : int 
            Maximum dipole offset in time bins to be allowed.
        lofar_centered : bool
            Beamforming station with respect to the LOFAR center (default) 
            or station centered.
        reftime : float
            Reference start time of the ob. computed for all stations.
        dm : float
            Dispersion measure used to freeze the data.
        overwrite : bool
            Remove output files if they previously exist.
        test : bool
            Only a few subbands are beamformed to check that there are 
            no errors.
        """
        
        # Input h5py data
        self.infile = infile #h5py.File(infile,'r')

        # Output filename
        self.bffilename = bffilename # TODO: Replace by bffile

        # Coordinates
        self.ra = ra
        self.dec = dec

        # Reference time
        if reftime:
            self.reftime = reftime
            self.time_key = int(np.floor(reftime))
            self.time_hr = datetime.fromtimestamp(reftime).strftime(
                    "%Y-%m-%d %H:%M:%S")
            self.sample_number = int((reftime%1)/(5e-9))
            self.slice = (reftime%1)

        # Other inputs
        self.overwrite = overwrite
        self.test = test
        self.pol = pol
        self.substation = substation
        self.offset_max_allowed = offset_max_allowed
        self.lofar_centered = lofar_centered
        self.dm = dm
        print("LOFAR centered:", lofar_centered)

    #-------------------------------------------------------------------
    # Internal functions called by the beamformer

    def __open_tbb_files(self):
        """Opening TBB files"""

        self.tbb_files = []
        for f in self.infile :
            try:
                tbb_file = h5py.File(f, 'r')
                self.tbb_files.append(tbb_file)
                print(tbb_file[tbb_file.keys()[0]])
            except:
                print("Unable to open file ", f)


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
            if k not in ['FILENAME', 'FILEDATE'] :
                attributes = []
                for f in self.tbb_files:
                    try:
                        attributes.append(f.attrs[k])
                    except:
                        print("Failed to open attribute ",k, "in file ",f)
                try:
                    attributes = []
                    if not all(attr == self.tbb_files[0].attrs[k] 
                           for attr in attributes):
                        print("ERROR: Attribute ",k, " not matching every file")
                        print([f.attrs[k] for f in self.tbb_files])
                        sys.exit()
                except:
                    if not all(attributes[i] == att for f in self.tbb_files 
                           for i,att in enumerate(self.tbb_files[0].attrs[k])):
                        print("ERROR: Attribute ",k, " not matching every file")
                        print([f.attrs[k] for f in self.tbb_files])
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
            caltabledir='/data/holography/Holog-20191212-1130/caltables/',
            parsetname='/home/veen/scripts/alert/StationCalibration.parset'):
        """Reading calibration tables"""

        # Reading calibration table
        filter_selection = (self.tbb_files[0].attrs[u'FILTER_SELECTION']
                .replace('A_','A-'))
        caltablename = (caltabledir + "CalTable-" + "" 
                + "{:03d}".format(self.station_number) 
                + '-' + filter_selection + '.dat')
        caltable = ct.readTable(caltablename)[1]
        
        # Reading clock offset
        self.clock_offset = float(md.getClockCorrectionParset(parsetname,
                self.station_name, antennaset=self.substation))
 
        return caltable        

    def __dipole_in_selection(self, dipole):
        """Checking if dipole belongs to substation"""

        d=int(dipole[-3:])
        if self.substation=='HBA0':
            if d%2 == self.pol and d<48: return True
        if self.substation=='HBA1':
            if d%2 == self.pol and d>=48: return True
        if self.substation=='HBA':
            if d%2 == self.pol: return True
        return False

    def __dipoles_subbands(self):
        """Function selecting dipoles and subbands from input data"""

        s = self.station
        # Selecting dipoles
        self.dipoles = [[str(d) for d in f[s].keys()] for f in self.tbb_files]
        self.selected_dipoles = [d for dip in self.dipoles for d in dip 
                if self.__dipole_in_selection(d)]
        self.selected_dipoles.sort()
        self.ndipoles = len(self.selected_dipoles)


        # Checking if there are dipoles
        if self.ndipoles == 0:
            print("ERROR: No dipoles selected for",self.station,self.substation)
            sys.exit() 

        # Select subbands
        self.subbands = np.unique([str(sb) for sb in f[s][d].keys() 
                for f in self.tbb_files for d in f[s].keys()])
        self.subbands.sort()
        self.nsubbands = len(self.subbands)

        # Define datalength
        d0 = self.tbb_files[0][s].keys()[0]
        sb0 = self.tbb_files[0][s][d0].keys()[0]
        #TODO: check if this value can we obtained from the header
        self.datalength = len(self.tbb_files[0][s][d0][sb0]) #585000
        print("datalength ", self.datalength)

        # Define bandwidth
        self.bandwidth = (self.tbb_files[0][s][d0][sb0].attrs['BANDWIDTH'])
        self.bandwidth_unit = (self.tbb_files[0][s][d0][sb0]
                .attrs['BANDWIDTH_UNIT'])

        #return dipoles, subbands

    def __subband_offsets(self, sb):
        """Computing time offsets for each dipole of the subband"""

        s = self.station

        # Select all dipoles with the given subband for the given polarisation
        available_dipoles = [str(d) for f in self.tbb_files for d in f[s].keys()
                if sb in f[s][d] 
                if d in self.selected_dipoles]
        available_dipoles.sort()
     
        # Calculating start time from offset
        d0 = available_dipoles[0]
        num = [i for i,dip in enumerate(self.dipoles) for d in dip if d0==d][0]

        frequency = self.tbb_files[num][s][d0][sb].attrs['CENTRAL_FREQUENCY']
        timeresolution = self.tbb_files[num][s][d0][sb].attrs['TIME_RESOLUTION']

        dm_delay = dd.freq_to_delay(self.dm, np.array([frequency,2e9]))[0]
        sb_starttime = self.reftime + self.clock_offset + dm_delay
 
        # Starttimes from files
        starttimes = [f[s][d][sb].attrs['TIME']+
                timeresolution*f[s][d][sb].attrs['SLICE_NUMBER']
                for i,f in enumerate(self.tbb_files)
                for d in self.dipoles[i] if d in available_dipoles]
        datalengths = [f[s][d][sb].shape[0]
                for i,f in enumerate(self.tbb_files)
                for d in self.dipoles[i] if d in available_dipoles]

        if min(starttimes) > sb_starttime:
            print("ERROR: min starttime > subband reftime", s, sb)
            #sys.exit()

        minstarttime = (int(np.floor(min(starttimes))), 
                int(min(starttimes)%1/timeresolution))
        maxstarttime = (int(np.floor(max(starttimes))), 
                int(max(starttimes)%1/timeresolution))
        diffstarttimes = (maxstarttime[0] - minstarttime[0] 
                + maxstarttime[1] - minstarttime[1]) / timeresolution

        # Calculating offsets
        offsets2  = [int(ceil((st - sb_starttime) / timeresolution))
                for st in starttimes]
        flag_offsets = [i for i,o in enumerate(offsets2) if o > 0]
        available_dipoles = [d for i,d in enumerate(available_dipoles)
                if i not in flag_offsets]
        offsets = [-o for i,o in enumerate(offsets2) if i not in flag_offsets]
        datalength = min(datalengths) - max(offsets)

        return offsets, available_dipoles, datalength

    def __subband_beamforming(self, 
            sb, offsets, available_dipoles, datalength, caltable):
        """
        Beamforming one subband

        Description
        -----------
        In the beamforming we correct for the calibration delay 
        and the beamforming weights.
        We're weighing by the number of available dipoles, 
        as this can be different per dipole

        Parameters
        ----------
        sb : 
            subband id.
        offsets :
            given by __subband_offsets
        available dipoles:
            dipoles with data for the given polarization
        datalength:
            number of time bins
        caltable:
            calibration table given by read_calibration

        Output
        ------
        Creation of subband dataset

        """

        s = self.station
        sbnr = int(sb[-3:])

        bfdata = np.zeros(shape = (datalength,), dtype = np.complex)
        sbdata_complex = np.zeros(shape = (datalength,), dtype=np.complex)
        d0 = available_dipoles[0]
        #num = np.where(np.char.find(self.dipoles, d0)!= -1)[0][0]
        num = [i for i,dip in enumerate(self.dipoles) for d in dip if d0==d][0]

        # Calculate time for the weights.
        # TODO: Is it subband time or reftime?
        t = utc2jd(self.tbb_files[num][s][d0][sb].attrs['TIME'])
        weights = azel2beamweights(getRADEC_2_AZEL(self.ra,self.dec,t), 
                self.station_name, 
                self.tbb_files[num][s][d0][sb].attrs[u'CENTRAL_FREQUENCY'],
                self.substation, lofarcentered=self.lofar_centered).toNumpy()
        ndipoles = len(available_dipoles) 

        # The beamforming part
        sbweight = np.zeros(shape = datalength)
        for d,o in zip(available_dipoles, offsets):
            #board = np.where(np.char.find(self.dipoles, d) != -1)[0][0]
            board = [i for i,dipgroup in enumerate(self.dipoles) 
                    for dip in dipgroup if d==dip][0]
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

        # Zero-padding
        pad = self.datalength-datalength
        bfdata = np.append(bfdata, np.zeros(pad))
        sbweight = np.append(sbweight, np.zeros(pad))

        # Create subband dataset
        self.bffile[s]['BFDATA'].create_dataset(sb,data=bfdata)
        self.bffile[s]['WEIGHTS'].create_dataset(sb,data=sbweight)


    def __subband_metadata(self, sb, offsets, available_dipoles):
        """Adding subband metadata"""

        s = self.station
        d0 = available_dipoles[0]
        #num = np.where(np.char.find(self.dipoles, d0)!= -1)[0][0]  
        num = [i for i,dip in enumerate(self.dipoles) for d in dip if d0==d][0]
        
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
            'FILEDATE': self.filedate,

            # Keywords from data attributes
            'TIME_RESOLUTION': fds.attrs['TIME_RESOLUTION'],
            'BANDWIDTH_UNIT': fds.attrs['BANDWIDTH_UNIT'],
            'SAMPLES_PER_FRAME': fds.attrs['SAMPLES_PER_FRAME'], 
            'TIME_RESOLUTION_UNIT': fds.attrs['TIME_RESOLUTION_UNIT'],
            'DATA_LENGTH': fds.attrs['DATA_LENGTH'], 
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
            'SUBSTATION': self.substation,
            'NOF_SELECTED_DATASETS': self.ndipoles,
            'NOF_DIPOLE_DATASETS': self.ndipoles,
            'DIPOLE_NAMES': 
                    [str(d) for d in self.selected_dipoles],
            'SELECTED_DIPOLES': 
                    [str(d.replace('DIPOLE_','')) for d
                    in self.selected_dipoles],
            'SELECTED_DIPOLES_INDEX': range(self.ndipoles),
            'CHANNEL_ID': 
                    [int(d.replace('DIPOLE_','')) for d 
                    in self.selected_dipoles],
            'SUBBANDS': [str(sb) for sb in self.subbands],
            'NOF_SUBBANDS': self.nsubbands,
            'SAMPLE_FREQUENCY': 200000000.0,
            'SAMPLE_FREQUENCY_UNIT': f.attrs["CLOCK_FREQUENCY_UNIT"], 
            'SAMPLE_FREQUENCY_VALUE': 200.0,
            'BANDWIDTH': self.bandwidth,
            'BANDWIDTH_UNIT': self.bandwidth_unit,
            'TIME': self.time_key,
            'TIME_HR': self.time_hr,
            'SAMPLE_NUMBER': self.sample_number, 
            'SLICE_NUMBER': int(self.slice/fds.attrs['TIME_RESOLUTION']),

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

    def beamforming(self):
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
        self.filedate = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

        # Opening input data
        self.__open_tbb_files()

        # Defining stations
        self.__station_selection() 
        s = self.station
        station_name = self.station_name
        station_number = self.station_number

        # Creating output data groups
        self.__data_groups()

        # Read calibration tables and clock offsets 
        caltable = self.read_calibration()

        # Selecting dipoles and subbands
        self.__dipoles_subbands()

        # Defining subbands to loop over. Less subbands for a test
        if self.test :
            sbs = self.subbands[0:2]
            self.subbands = self.subbands[0:2]
            self.nsubbands = len(self.subbands)
        else :
            sbs = self.subbands[0:]

        # Loop over subbands
        for sb in sbs:

            print ("{} POL{} {} Subband {}/{}".format(self.station, self.pol,
                    self.substation, sb, sbs[-1]))
            sys.stdout.write("\033[F") # Cursor up one line

            # Select all dipoles with the given subband and given polarisation
            available_dipoles = [d for f in self.tbb_files for d in f[s].keys()
                    if sb in f[s][d].keys() if d in self.selected_dipoles]
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
        antenna_set = self.substation
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

            'antenna_set': antenna_set,
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

        self.__bf_dictionary()

        # Writing data        
        print("Writing FFT", self.station, "with shape", self.fftdata.shape)
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
            TIME = [f.attrs["TIME"]] * ndipoles,
            TIME_HR = [f.attrs["TIME_HR"]] * ndipoles,
            SAMPLE_NUMBER = [f.attrs["SAMPLE_NUMBER"]] * ndipoles,
            SLICE_NUMBER = [f.attrs["SLICE_NUMBER"]] * ndipoles,
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
            clock_offset = [float(md.getClockCorrectionParset(
                    '/home/veen/scripts/alert/StationCalibration.parset',
                    self.station_name, self.substation))] * ndipoles,
            MAXIMUM_READ_LENGTH = 204800000, # Check
            FFTSIZE = self.fftdata.shape[-1],
            NYQUIST_ZONE = [2] * ndipoles,
            FREQUENCY_RANGE = "", # [(100000000.0, 200000000.0)] * 96
            PIPELINE_VERSION = "UNDEFINED",

            # Derived keywords that are less relevant
            TIMESERIES_DATA = "", # Check
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
        
        Description
        -----------
        A .beam pycrtools file with the dynamic spectrum is created from the 
        .h5 TBB beamformed data.

        Parameters
        ----------
        dynspecfile : str
            dynamic spectrum output file
        DM : int
            dispersion measure of the source
        nch : int
            number of channels
        t_int : int 
            integration time
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
        self.substation = self.bffile.attrs["SUBSTATION"]
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
# Beamforming across stations
#------------------------------------------------------------------------

def add_stations(filenames, outbfdir, tbin=12, dm=0, incoherent=True, 
        skip_stations=[]):
    """
    Using pycrtools function addBeams to beamform across stations.
    
    Parameters
    ----------
    files : list 
        list of .beam station files
    tbin : int 
        bin time
    dm : float
        dispersion measure (0 if the data is already dedispersed)
    incoherent : bool
        True if incoherent beam addition, False if coherent 
    skip_stations : list 
        list of stations to skip when beamforming across stations
    """
    
    # Opening files
    #files_all = glob.glob(filenames)
    #files_all.sort()
    files = []

    for f in filenames:
      if len(skip_stations) > 0:
        st = f.split('/')[-1].split('_')[1]
        if st in skip_stations:
          print("Skipping station", st)
          continue
      try:
        cr.open(f)
        files.append(f)
      except:
        print(f, "could not be opened")

    files.sort()
    print("Reading files", files)

    # Beamforming
    beams = cr.open(files)

    TAB, dynspec, cleandynspec = bt.addBeams(beams, dyncalc=True, tbin=tbin,
            dm=dm, clean=True, incoherent=incoherent)

    # Defining output names
    f = files[0].split('/')[-1].replace('.beam', '').split('_')
    obs = [l for l in f if 'L' in l][0]
    pol = [p for p in f if 'pol' in p][0]
    #hba = [h for h in f if 'HBA' in h][0]
    date = [d for d in f if 'D' in d and 'T' in d][0]

    if outbfdir[-1] != '/': outbfdir = outbfdir + '/'

    tabname = outbfdir+'{0}_{1}_TAB_{2}'.format(obs, date, pol)
    dsname  = outbfdir+'{0}_{1}_dynspec_{2}'.format(obs, date, pol)
    cdsname = outbfdir+'{0}_{1}_cleandynspec_{2}'.format(obs, date, pol)

    print("Writing files:")
    print(tabname, '\n', dsname, '\n', cdsname)

    # Saving beamformed data
    TAB.write(tabname + '.beam')
    dynspec.write(dsname + '.beam')
    cleandynspec.write(cdsname + '.beam')

    npTAB = TAB.toNumpy()
    npdynspec = dynspec.toNumpy()
    npcleandynspec = cleandynspec.toNumpy()

    np.save(tabname, npTAB)
    np.save(dsname,  npdynspec)
    np.save(cdsname, npcleandynspec)

if __name__=='__main__':

    # Command line options
    parser = OptionParser()

    parser.add_option("-f", "--outfiledir", type="str",
           default="/data/projects/COM_ALERT/pipeline/analysis/marazuela/data/",
           help="Directory where the output files will be generated.",
           dest="outfiledir")
    parser.add_option("--dm", "--dispersion_measure", dest="dm", type="float",
           default=0, 
           help="Dispersion measure difference between target and freezing DM.")
    parser.add_option("-b", "--tbin", "--time_bin", dest="tbin",
           type="int", default=12,
           help="Maximum offset between subbands in time bins.")
    parser.add_option("-i", "--incoherent", dest="incoherent",
           help=("If provided, incoherent beamforming. Default: Coherent"),
           action="store_true", default=False)
    parser.add_option("-s", "--skip_stations", type="str",
           default="", help="Stations to skip when adding them together",
           dest="skip_stations")

    (options, args) = parser.parse_args()

    if options.outfiledir[-1] != '/': options.outfiledir = options.outfiledir + '/'

    filenames = []
    for f in args[:]: filenames.append(glob.glob(f)[0])
    filenames.sort()

    skip_stations = [st[:5] for st in options.skip_stations.upper().split(',')]

    # Beamforming across stations
    add_stations(filenames, options.outfiledir, tbin=options.tbin, 
           dm=options.dm, incoherent=options.incoherent, 
           skip_stations=skip_stations)

#------------------------------------------------------------------------
# End of the script
