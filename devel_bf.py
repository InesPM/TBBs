from __future__ import print_function
from __future__ import division
import h5py
import numpy as np
import caltable as ct
from radec2azel import *
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
    
    def __init__(self, infile, outfile):
        self.infile = infile
        self.outfile = outfile

    def beamforming(self):
        #----------------------------------------------------------------
        # Opening files

        # Open given tbb file
        f=h5py.File(self.infile,'r')

        # Open output file
        f_out=h5py.File(self.outfile,'w')

        #----------------------------------------------------------------
        # Select stations in file. We now only select the first station. 
        # We could loop over stations, but at the moment one station is written per file
        stations=[k for k in f.keys() if 'STATION' in k]
        station=stations[0]
        station_name=station[-5:]
        # abbreviation for station, TODO: expand in final code
        s=station
        # station number
        st_nr=int(s[-3:])
        

        # Creating keywords output file
        print("Creating hdf5 output file")
        for k in f.attrs.keys():
            try:
                f_out.attrs.create(k,data=f.attrs[k])
            except:
                if k in ["TARGETS","OBSERVATION_STATIONS_LIST"]:
                    f_out.attrs.create(k,data=[vv for vv in f.attrs[k]])
                else:
                    print("failed to create key",k)
        f_out.create_group(s)
        f_out[s].create_group('WEIGHTS')
        f_out[s].create_group('BFDATA')
        for k in f[s].attrs.keys():
            try:
                f_out[s].attrs.create(k,data=f[s].attrs[k])
            except:
                print("failed to create key",k)
        # Additional keywords
        for k in ["DIPOLE_NAMES","SELECTED_DIPOLES"]:
            try :
                f_out.attrs.create(k,data=[str(d.replace('DIPOLE_', '')) for d in f[station].keys() if int(d[-3:])%2==pol])
            except :
                print("failed to create key",k)

        # Read calibration tables. These should also be in a directory in the image, but can be old.
        filter_selection=f.attrs[u'FILTER_SELECTION'].replace('A_','A-')
        caltabledir='/data/holography/Holog-20191212-1130/caltables/'
        caltablename=caltabledir+"CalTable-"+""+"{:03d}".format(st_nr)+'-'+filter_selection+'.dat'
        caltable=ct.readTable(caltablename)[1]
        print("Stations present",stations," selected",station)

        # Selecting dipoles
        print("Selection dipoles")
        dipoles=f[station].keys()

        # Select subbands
        subbands=np.unique(np.concatenate([f[station][d].keys() for d in dipoles]))
        subbands.sort()

        # TODO loop over subbands, first a porton, then all
        subband=subbands[0]
        for sb in subbands[0:5]:
            sbnr=int(sb[-3:])
            # Select all dipoles that have this subband in their keys, for the even or odd polarisations
            available_dipoles=[d for d in dipoles if sb in f[s][d].keys() if int(d[-3:])%2==pol ]
            available_dipoles.sort()
            if len(available_dipoles)==0:
                continue


            print("Calculating offsets")
            # use one dipole to select some of the keywords. We assume they are even over the subband
            d=available_dipoles[0]
            timeresolution=f[s][d][sb].attrs['TIME_RESOLUTION']
            starttimes=[(f[s][d][sb].attrs['TIME'],timeresolution*f[s][d][sb].attrs['SLICE_NUMBER']) for d in available_dipoles]
            datalengths=[f[s][d][sb].shape[0] for d in available_dipoles]
            #assert min(datalengths)==max(datalengths) # time to develop a flagging function
            minstarttime=min(starttimes)
            maxstarttime=max(starttimes)
            diffstarttimes=(maxstarttime[0]-minstarttime[0]+maxstarttime[1]-minstarttime[1])/timeresolution
            offsets2=[int(math.ceil(((st[0]-minstarttime[0])+(st[1]-minstarttime[1]))/timeresolution)) for st in starttimes]
            flag_offsets=[num for num,o in enumerate(offsets2) if o>offset_max_allowed]
            available_dipoles=[d for num,d in enumerate(available_dipoles) if num not in flag_offsets]
            starttimes=[(f[s][d][sb].attrs['TIME'],timeresolution*f[s][d][sb].attrs['SLICE_NUMBER']) for d in available_dipoles]
            minstarttime=min(starttimes)
            maxstarttime=max(starttimes)
            offsets=[int(round((maxstarttime[0]-st[0]+maxstarttime[1]-st[1])/timeresolution)) for st in starttimes]
            datalength=min(datalengths)-max(offsets)

            print("subband",sb,"# Available dipoles",len(available_dipoles))

            print("Creating output datasets")
            bfdata=np.zeros(shape=(datalength,),dtype=np.complex)
            sbdata_complex=np.zeros(shape=(datalength,),dtype=np.complex)
            d0=available_dipoles[0]
            # Calculate time for the weights. We could vary the weight over the observation if this is longer than one second
            t=utc2jd(f[s][d0][sb].attrs['TIME'])
            weights=azel2beamweights(getRADEC_2_AZEL(ra,dec,t),station_name,f[s][d0][sb].attrs[u'CENTRAL_FREQUENCY']).toNumpy()
            ndipoles=len(available_dipoles) 
            # The beamforming part
            sbweight=np.zeros(shape=datalength)
            for d,o in zip(available_dipoles,offsets):
                if d in ['DIPOLE_145000000','DIPOLE_145000002']:
                    continue
                print("Analysing dipole",d,"at offset",o,)
                dnr=int(d[-3:])
                # Read data from the offset, so that all are equal
                sbdata=f[s][d][sb][o:o+datalength]
                # Convert to complex numbers
                try:
                    sbdata_complex.real=sbdata['real']
                    sbdata_complex.imag=sbdata['imag']
                except:
                    sbdata_complex=sbdata
                # In the beamforming we correct for the calibration delay and the beamforming weights
                # We're weighing by the number of available dipoles, as this can be different per dipole
                # TODO: Check signs
                dipoledata=sbdata_complex*caltable[dnr,sbnr]*weights[dnr]
                dipoleweight=dipoledata/dipoledata
                dipoleweight[np.isnan(dipoleweight)]=0
                sbweight+=dipoleweight
                nzeroes=np.sum(np.abs(dipoledata)==0)
                bfdata+=dipoledata
                print("nzeroes",nzeroes/len(bfdata))

            sbweight[sbweight==0]=1
            bfdata/=np.sqrt(sbweight)
            # Create subband dataset
            f_out[s]['BFDATA'].create_dataset(sb,data=bfdata)
            f_out[s]['WEIGHTS'].create_dataset(sb,data=sbweight)

      
            # Add metadata
            for k in list(f[s][d][sb].attrs):
                if k not in ['FLAG_OFFSETS','SLICE_NUMBER']:
                    f_out[s]['BFDATA'][sb].attrs.create(k,data=f[s][d][sb].attrs[k])
            # Correct slice number for the offsets
            # TODO: Add check that all offsets are the same
            f_out[s]['BFDATA'][sb].attrs.create('SLICE_NUMBER',data=f[s][d][sb].attrs['SLICE_NUMBER']+np.array(offsets)[np.array(available_dipoles)==d][0])
            for k in [
              #u'TILE_BEAM', # TODO- read all TILE keys
              #u'TILE_BEAM_UNIT',
              u'STATION_ID',
              #u'TILE_BEAM_FRAME',
              #u'SAMPLE_FREQUENCY_UNIT',
              u'NYQUIST_ZONE',
              u'SAMPLE_FREQUENCY']:
                f_out[s]['BFDATA'][sb].attrs.create(k,data=f[s][d].attrs[k])
            for k in [
              u'RCU_ID',
              u'RSP_ID']:
                f_out[s]['BFDATA'][sb].attrs.create(k,data=[f[s][d].attrs[k] for d in available_dipoles])
            for k in [u'DIPOLE_IDS']:
                f_out[s]['BFDATA'][sb].attrs.create(k,data=[str(d) for d in available_dipoles])

        f_out.close()

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
offset_max_allowed = options.offset_max_allowed
outfilename = options.outfiledir + ffull.split('/')[-1][0:-3]+'bf_pol'+str(pol)+'.h5'

assert(pol in [0,1])

#------------------------------------------------------------------------
# Beamforming
#------------------------------------------------------------------------

b = beamformer(ffull, outfilename)
b.beamforming()

print("File written: ",outfilename)
