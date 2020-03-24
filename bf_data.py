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
import tbb_beamformer as bf

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
assert(pol in [0,1])

offset_max_allowed = options.offset_max_allowed

# Defining output file names
bffilename, dsfilename = bf.outfilenames(options.outfiledir, files, pol)

#------------------------------------------------------------------------
# Beamforming
#------------------------------------------------------------------------

b = bf.BeamFormer(infile=files, bffilename=bffilename, pol=pol, 
    offset_max_allowed=offset_max_allowed)
b.beamforming()
#b.convert2beam(dsfilename)


print("File written: ", bffilename)
print("File written: ", dsfilename)
