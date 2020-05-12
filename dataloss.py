from __future__ import print_function
from __future__ import division
import sys
import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

def dataloss(fn):
    #path = '/data/projects/COM_ALERT/pipeline/analysis/marazuela/data/'
    #out = path + fn.split('/')[-1].replace('.h5', '.png')
    #title = fn.split('/')[-1].replace('.h5', '')

    ff = h5py.File(fn)

    station = str(ff.keys()[0])
    dipoles = ff[station].keys()

    nzeros  = 0
    datalen = 0

    #data_arr = [ff[station][kk][jj][:].real for kk in ff[station].keys() 
    #        for jj in ff[station][kk].keys()]
    #data_arr = np.concatenate(data_arr)

    zeros = np.array((0, 0), dtype=[('real', '<i2'), ('imag', '<i2')])
    #np.load('/home/marazuela/scripts/zeros.npy')

    for d in dipoles:
        subbands = ff[station][d].keys()
        for sb in subbands:
            datalen += len(ff[station][d][sb][:])
            nzeros  += len(np.where(ff[station][dipoles[0]][subbands[0]][:]==zeros)[0])

    dataloss = round(100.0 * nzeros / datalen, 1)

    return dataloss

if __name__=='__main__':
    
    # Command line arguments
    parser = OptionParser()
    (options, args) = parser.parse_args()

    # Input files
    files = []
    stations = []
    for a in args[:]:
        f = glob.glob(a)[0]
        files.append(f)
        stations.append(f.split('/')[-1].split('_')[1])

    files.sort()
    stations = np.unique(stations)
    stations.sort()

    for s in stations:
        print('\n', s, end='\t')
        for f in files:
            if s in f:
                dl = dataloss(f)
                print(dl, end='\t')           






