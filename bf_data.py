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
import os
import glob
from optparse import OptionParser
import tbb_beamformer_devel as bf

#------------------------------------------------------------------------
# Defining functions
#------------------------------------------------------------------------

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

def outfilenames(outfiledir, files, pol, subs, test=False):
    """Defining output file names"""

    # Checking number of input files
    if len(files) > 1:
        fname = findstem(files).split('/')[-1] + 'tbb'
    else :
        fname = files[0].split('/')[-1][0:-3]

    # Defining f=general output directory
    stem = fname.split('_')    
    obsdir = stem[0] + ('_') + stem[2].split('.')[0]
    if test:
        obsdir = obsdir + '_test'

    # Creating output directories
    try:
        os.mkdir(outfiledir + obsdir)
        print("Directory", outfiledir + obsdir, "created")
    except:
        print("Directory", outfiledir + obsdir, "already exists")

    bffiledir = outfiledir + obsdir + '/bfh5_data/'
    try:
        os.mkdir(bffiledir)
        print("Directory", bffiledir, "created")
    except:
        print("Directory", bffiledir, "already exists")

    beamfiledir = outfiledir + obsdir + '/beam_data/'
    try:
        os.mkdir(beamfiledir)
        print("Directory", beamfiledir, "created")
    except:
        print("Directory", beamfiledir, "already exists")

    # Defining file names
    outbf = bffiledir + fname + '_bf_'+str(subs)+'_pol'+str(pol)+'.h5'
    outdynspec = beamfiledir + fname + '_fft_'+str(subs)+'_pol'+str(pol)+'.beam'
    return outbf, outdynspec

#------------------------------------------------------------------------
# Command line options
#------------------------------------------------------------------------
parser = OptionParser()

parser.add_option("-r", "--right_ascension", dest="ra", type="float", 
        default=0.92934186635, 
        help="Right ascension of the source in degrees. Default RA: B0329")
parser.add_option("-d", "--declination", dest="dec", type="float", 
        default=0.952579228492, 
        help="Declination of the source. Default Dec: B0329")
parser.add_option("-p", "--polarisation", dest="pol", type="int", default=0,
        help="Polarisation to beamform.")
parser.add_option("-s", "--substation", dest="substation", type="str", 
        default='HBA0', help="Polarisation to beamform.")
parser.add_option("-o", "--offset_max_allowed", dest="offset_max_allowed", 
        type="int", default=400, 
        help="Maximum offset between subbands in time bins.")
parser.add_option("-f", "--outfiledir", type="str", 
        default="/data/projects/COM_ALERT/pipeline/analysis/marazuela/data/",
        help="Directory where the output files will be generated.", 
        dest="outfiledir")
parser.add_option("-t", "--test", action="store_true", 
        help="Test script by beamforming with 2 subbands", dest="test")

(options, args) = parser.parse_args()

# Input file
files = []
for f in args[:]:
    files.append(glob.glob(f)[0])
#files=glob.glob(args[:])
files.sort()
files.insert(0, files.pop(-1))
print('Opening following files:')
for f in files : print(f)

# Checking if file exists or is a list
assert len(files) >= 1, "Input file(s) not found"

# Other inputs
(ra,dec) = (options.ra, options.dec)
pol = options.pol
subs = options.substation
assert(pol in [0,1])

offset_max_allowed = options.offset_max_allowed

# Defining output file names
bffilename, dsfilename = outfilenames(options.outfiledir, files, pol, subs, test=options.test)

#------------------------------------------------------------------------
# Beamforming
#------------------------------------------------------------------------

b = bf.BeamFormer(infile=files, bffilename=bffilename, pol=pol, 
    substation=options.substation, offset_max_allowed=offset_max_allowed, 
    test=options.test, overwrite=True)
b.beamforming()
b.convert2beam(dsfilename)

print("File written: ", bffilename)
print("File written: ", dsfilename)
