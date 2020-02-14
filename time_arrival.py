import sys
import itertools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors, colorbar

sys.path.insert(0, '/home/marazuela/scripts/tbb-tools')
import read_tbb_data

def dedisp(t, sb, dm=0) :
    "returns dedispersed time"
    f = []
    t_dd = []
    for i in range(len(sb)) :
        f.append([])
        t_dd.append([])
        for j in range(len(sb[i])) :
          nu = (1e8 + 195312.5 * int(sb[i][j])) * 1e-9  
          f[i].append(nu)
          td = t[i][j] - 4.15e-3 * dm * (1/nu**2 - 1/0.2**2)
          t_dd[i].append(td)
    return f, t_dd 

def plot_time_arrival(data, dm=0):

    # Getting data
    tarr, sbb = data.get_time_arr_dip(fn)
    t0 = min([t for times in tarr for t in times])

    title = fn.split('/')[-1].replace('.h5', '')
    path = '/home/marazuela/plots/'
    out = path + 'tarr_' + title + '.pdf'
    print(out)

    xl = (max(max(tarr)) - min(min(tarr)))

    c1 = cm.magma_r(np.linspace(0, 1, nsb))
    c2 = cm.viridis(np.linspace(0, 1, len(tarr)))


    # Plotting
    plt.title(title)
    plt.xlabel('Arrival time')
    plt.ylabel('Number')

    plt.xlim(min(min(tarr))-xl*0.1, max(max(tarr))+xl*0.1)
    plt.ylim(0, nsb)

    ####
    # Figure 1 : arrival time vs. subband

    plt.subplot(2,1,2)
    for i in range(len(tarr)) :
      f, t_dd = dedisp(tarr, sbb, dm)
      plt.scatter(t_dd[i], sbb[i], marker='o', c=c2[i], cmap='viridis', alpha=0.5, s=50)

    plt.xlim(min(min(t_dd))-xl*0.1, max(max(t_dd))+xl*0.1)

    plt.xlabel('Time (s)')
    plt.ylabel('Subband')

    plt.savefig(out, pad_inches=0, bbox_inches='tight')
    plt.show()

def plot_starting_0(data) :

    # TODO: add a data shape checker. For '/data/projects/COM_ALERT/tbb/L597863_CS007_D20200211T161560.000Z_R005_tbb.h5', station CS007 is doubled.

    title = fn.split('/')[-1].replace('.h5', '')
    path = '/home/marazuela/plots/'
    out = path + 't0_dipoles_' + title + '.pdf'
    print(out)

    # Getting data
    tarr, sbb = data.get_time_arr_dip(fn)
    t0 = min([t for times in tarr[0] for t in times])

    # Sorting subbands
    subbands = np.sort(np.unique([val for sublist in sbb[0] for val in sublist]))

    tarr_0 = np.zeros(len(subbands))

    for i in range(len(tarr[0])) :
        for j in range(len(tarr[0][i])) :
            index = np.where(subbands == sbb[0][i][j])
            if tarr[0][i][j] - t0 <= 40*5.12e-6 :
                tarr_0[index] += 1

    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    ax.plot(subbands, tarr_0)

    ax.set_ylim(0, max(tarr_0)*1.1)
    ax.set_xlabel("Subband number")
    ax.set_ylabel("Dipoles starting at 0")

    plt.savefig(out, pad_inches=0, bbox_inches='tight')
    #plt.show()


if __name__=='__main__':

    fn = sys.argv[1]
    #fn = '/data/projects/COM_ALERT/tbb/L597863_RS305_D20200128T153519.000Z_tbb.h5'
    nsb = int(sys.argv[2])
    dm = float(sys.argv[3])

    T = read_tbb_data.TBBh5_Reader(fn)
    #f, attr = T.read_h5()
    tarr, sbb = T.get_time_arr_dip(fn)
    t0 = min([t for times in tarr for t in times])    

    #t = []
    #for i in range(nsb) :
    #  t.append([])

    #tmin = np.amin(tarr)
    #for d in range(len(tarr)) :
    #  for i in range(nsb) :
    #    for j in range(len(sbb[d])) :
    #      if i+(400-nsb) == sbb[d][j] :
    #        t[i].append(tarr[d][j])
    #t = np.array(t)

    #xl = (max(max(tarr)) - min(min(tarr)))
 
    #c1 = cm.magma_r(np.linspace(0, 1, nsb))
    #c2 = cm.viridis(np.linspace(0, 1, len(tarr)))


    ####
    # Figure 2 : histogram 
      
    #plt.subplot(2,1,1)

    #for i in range(nsb) :
    #  plt.hist(t[i], bins=(min(min(tarr)), max(max(tarr)), 20), range=(0, xl), edgecolor=c1[i], fc='None', linestyle='-', histtype='stepfilled', lw=2)

    #plt.yscale('log')    
    #plot_time_arrival(T)
    plot_starting_0(T) 



