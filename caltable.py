import struct
import numpy as np
import os,sys
import matplotlib.pyplot as plt


def readTable(filename):
    """Read phase calibration data for a station.

    *filename* filename of the caltable

    returns weights for 512 subbands.

    Examples (also for doctests):

    >>> pycr_metadata.station_phase_calibration_factors("CS302","LBA_OUTER")

    >>> pycr_metadata.station_phase_calibration_factors(122,"LBA_OUTER")

    """

    file = open(filename, 'rb')
    filesize = os.path.getsize(filename)
    if filesize < 2 * 786432:
        datasize = 786432
        nant = 96
        # reading in 96 * 512 * 2 doubles
        fmt = '98304d'
    else:
        datasize = 786432 * 2
        nant = 192
        # international stations 192 * 512 * 2 doubles
        fmt = '196608d'

    # We skip the header in the datafile
    # file.seek(filesize-datasize)

    #header = file.read(filesize - datasize)

    # Calculate the size to be read
    sz = struct.calcsize(fmt)

    # read from the file
    rawdata = file.read()
    header_end=str.find(rawdata,'HeaderStop\n',0,len(rawdata))+len('HeaderStop\n')
    header=rawdata[:header_end]
    # unpack so the data is represented as doubles
    data = struct.unpack(fmt, rawdata[header_end:header_end+sz])

    #
    data = np.array(data)
    data.resize(512, nant, 2)

    complexdata = np.empty(shape=(512, nant), dtype=complex)
    complexdata.real = data[:, :, 0]
    complexdata.imag = data[:, :, 1]

    # swap frequency and time axis
    return header, complexdata.transpose()


def writeTable(header, table, filename, overwrite=False):
    if not overwrite:
        assert not os.path.isfile(filename)
    file=open(filename,'w')
    if not checkHeader(header):
        raise ValueError,"Headers should start with Headerstart and end with HeaderStop\\n. Are you sure you are reading in the correct header and data?"
    nant=table.shape[0]
    if nant not in [96,192]:
        assert False # wrong number of antennas
    # swap frequency and time axis
    complexdata=table.transpose()

    # Create float array
    data=np.zeros(shape=(512,nant,2),dtype=np.float64)
    data[:,:,0]=complexdata.real
    data[:,:,1]=complexdata.imag
    data=list(data.flatten())

    # pack data and write to binary file
    if nant==96:
        fmt='98304d'
    elif nant==192:
        fmt = '196608d'
    sz=struct.calcsize(fmt)
    rawdata=struct.pack(fmt,*data)
    file.write(header)
    file.write(rawdata)
    file.close()


def plotTable(table,antlist=np.arange(0,192)):
    phaseticklabels = np.arange(-1.0, 1.5, 0.5)  # added by MHD
    phasetickpositions = phaseticklabels*np.pi  # added by MHD
    phaseticklabels = ['%.2f'%(x)+'$\pi$' for x in phaseticklabels]  # added by MHD
    for ant in antlist:
        plt.subplot(2,1,1)
        plt.plot(np.angle(table[ant]))
        plt.ylim(-np.pi, np.pi)  # added by MHD: restrict plot axes to -pi, +pi
        plt.yticks(phasetickpositions, phaseticklabels)  # added by MHD
        plt.subplot(2,1,2)
        plt.plot(np.abs(table[ant]))
    plt.xlabel('subband')  # MHD: changed from 'antenna' to 'subband'
    plt.ylabel('amplitude')
    plt.subplot(2,1,1)
    plt.ylabel('phase')

def calcDelay(table,antlist=np.arange(0,192)):
    # Add check to see if deltaphi is linear ?
    delays=np.median(np.diff(np.angle(table[:]), axis=1) / 195312.5 / 2 / np.pi, axis=1)
    return delays[antlist]

def addDelay(table,ant,delay):
    table=np.copy(table)
    subbandwidth=195312.5 # Hz
    nsubbands=512
    dFreq=np.arange(0,subbandwidth*nsubbands,subbandwidth)
    phaseoffsets=dFreq*2*np.pi*delay*1j
    phases=np.exp(phaseoffsets)
    table[ant]*=phases
    return table

def plotDelay(delays,antlist=np.arange(0,192)):
    plt.plot(antlist,delays)
    plt.title('Delays')
    plt.xlabel('Antenna ID')
    plt.ylabel('Delay [s]')

def checkHeader(header):
    if header[0:11]!='HeaderStart': return False
    if header[-11:]!='HeaderStop\n': return False
    return True




