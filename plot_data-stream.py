import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

fn = sys.argv[1]
station = sys.argv[2]
title = sys.argv[3]

path = '/data/projects/COM_ALERT/pipeline/analysis/marazuela/data/'
out = path + fn.split('/')[-1].replace('.h5', '.png')

#ff=h5py.File('/data/projects/COM_ALERT/tbb/L597863_RS106_D20191023T095157.000Z_tbb.h5')

station='STATION_' + station

ff = h5py.File(fn)

data_arr = []

for kk in ff[station].keys():
    for jj in ff[station][kk].keys():
        data_arr.append((ff[station][kk][jj][:].real))
data_arr = np.concatenate(data_arr)
print "# zeros",np.sum((data_arr['real']==0)&(data_arr['imag']==0)),"out of",len(data_arr)," dataloss ",round(100.0*np.sum((data_arr['real']==0)&(data_arr['imag']==0))/len(data_arr),1),"%"


# Plots

#plt.figure(figsize=(21, 8))
#plt.plot(data_arr, color='teal')
#plt.xlabel('Time sample')
#plt.savefig(out, pad_inches=0, bbox_inches='tight')
#plt.show()

