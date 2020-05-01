#! /usr/bin/env python
"""
Getting reference time to align stations to be beamformed
"""

from optparse import OptionParser
import pycrtools as cr
from pycrtools import metadata as md
import h5py
import os; import sys; import glob; import time
import numpy as np
import dedispersion as dd

#########################################################################
# Defining functions

def getRef_time(filelist, outfiledir, DM=26.76,
        fullparsetname='/home/veen/scripts/alert/StationCalibration.parset',
        substation='HBA0', offset_max_allowed=400, verbose=True):
    """
    Gets the reference time for a list of tbb files 
    (station which start observing last).
    """

    start = time.time()
    reftime = np.zeros(len(filelist), dtype=np.float128)

    for i,file in enumerate(filelist):
        if 'CS032' in file:
            # Skipping this file for now
            print "Skipping file", file
            continue
        try:  
            print "Opening file", file
            tbb = h5py.File(file, 'r')
            pass
        except:
            print "Unable to open file", file
            continue

        station = tbb.keys()[0]
        station_name = station[-5:]
        s = station

        if 'CS' in station_name:
            substation = 'HBA0'
        elif 'RS' in station_name:
            substation = 'HBA'

        # Getting clock offset
        print "Getting clock offset", station
        if fullparsetname:
            f_clock_offset = float(md.getClockCorrectionParset(fullparsetname,
                    station_name, antennaset=substation))
        else:
            f_clock_offset = float(md.getClockCorrection(station, 
                    antennaset=substation))

        # Selecting first subband from first dipole
        dipoles = tbb[station].keys()
        subbands = np.unique([str(sb) for d in tbb[s].keys()
                for sb in tbb[s][d].keys()])
        timeres = tbb[s][dipoles[0]][subbands[0]].attrs['TIME_RESOLUTION']

        # Getting start times from file
        print "Getting start time", station
        starttimes = np.zeros((len(subbands),len(dipoles)),dtype=np.float128)
        frequency = np.zeros(len(subbands))
        sb2j = dict()

        for j,sb in enumerate(subbands):
          k = 0
          sb2j[sb]=j
          while sb not in tbb[s][dipoles[k]].keys(): k += 1
          frequency[j] = tbb[s][dipoles[k]][sb].attrs['CENTRAL_FREQUENCY']

        for k,d in enumerate(dipoles):
            for sb in tbb[s][d].keys():
              j = sb2j[sb]
              starttimes[j,k] = (tbb[s][d][sb].attrs['TIME'] 
                      + timeres*tbb[s][d][sb].attrs['SLICE_NUMBER'])

        starttimes[starttimes == 0] = np.nan

        # Time correction
        print "Time correction", station
        dm_delay = np.delete(dd.freq_to_delay(DM,np.append(frequency,2e9)),-1)
        dedisptime = (starttimes.T-dm_delay).T
        offsettime = dedisptime + f_clock_offset

        # Min, max time
        minstarttime = np.nanmin(offsettime)
        maxstarttime = np.nanmax(offsettime)
        minsbstarttime = np.nanmin(offsettime, axis=1)
        medsbstarttime = np.median(minsbstarttime)

        print "Reftime", station
        mask = np.where(offsettime>medsbstarttime+timeres*offset_max_allowed)
        if len(mask) > 0:
            reftime[i] = medsbstarttime + timeres * offset_max_allowed
        elif medsbstarttime + timeres * offset_max_allowed > maxstarttime:
            reftime[i] = maxstarttime
        else:
            print "Couldn't define reftime", tbb
        
        print "Elapsed time", time.time() - start, "s"

    # Ending loop
    if verbose:
        print 'Reference time: ' , np.max(reftime) , 'from file ' 
    # Saving data
    try:
        np.save(outfiledir+'reftime', max(reftime))
        np.save(outfiledir+'reftimes', reftime)
    except:
        print("Could not save reference times")

    return max(reftime), reftime

#########################################################################
# Main script

if __name__=='__main__':

    # Command line options
    parser = OptionParser()

    parser.add_option("-t", "--test", action="store_true", 
           help="Test script by beamforming with 2 subbands", dest="test")
    parser.add_option("--dm", "--dispersion-measure", dest="dm", type="float", 
           default=26.8, help="Dispersion measure.")
    parser.add_option("-o", "--offset_max_allowed", dest="offset_max_allowed", 
           type="int", default=2000, 
           help="Maximum offset between subbands in time bins.")
    parser.add_option("-f", "--outfiledir", type="str",
           default="/data/projects/COM_ALERT/pipeline/analysis/marazuela/data/",
           help="Directory where the output files will be generated.",
           dest="outfiledir")

    (options, args) = parser.parse_args()

    # Input file
    files = []
    for f in args[:]:
        files.append(glob.glob(f)[0])
    files.sort()

    if options.outfiledir[-1] != '/':
        options.outfiledir = options.outfiledir + '/'

    # Computing reference times
    reftime, reftimes = getRef_time(files, options.outfiledir, 
            DM=options.dm, offset_max_allowed=options.offset_max_allowed)

# End of script
