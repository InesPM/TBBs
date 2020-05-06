import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

fn = sys.argv[1]
#station = sys.argv[2]
#title = sys.argv[3]

path = '/data/projects/COM_ALERT/pipeline/analysis/marazuela/data/'
out = path + fn.split('/')[-1].replace('.h5', '.png')
title = fn.split('/')[-1].replace('.h5', '')

ff = h5py.File(fn)

station = str(ff.keys()[0])

data_arr = [ff[station][kk][jj][:].real 
           for kk in ff[station].keys() 
           for jj in ff[station][kk].keys()]
data_arr = np.concatenate(data_arr)

zeros = np.load('/home/marazuela/scripts/zeros.npy')
#print "# zeros",np.sum((data_arr['real']==0)&(data_arr['imag']==0)),"out of",len(data_arr)," dataloss ",round(100.0*np.sum((data_arr['real']==0)&(data_arr['imag']==0))/len(data_arr),1),"%"
#print "# zeros",np.sum((data_arr==zeros)),"out of",len(data_arr)," dataloss ",round(100.0*np.sum((data_arr==zeros))/len(data_arr),1),"%"
print round(100.0*np.sum((data_arr==zeros))/len(data_arr),1)

# Plots

#plt.figure(figsize=(21, 8))
#plt.plot(data_arr, color='teal')
#plt.xlabel('Time sample')
#plt.savefig(out, pad_inches=0, bbox_inches='tight')
#plt.show()

