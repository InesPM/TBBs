
import h5py
import numpy as np

# Input and output files.

input_file='/data/projects/COM_ALERT/pipeline/data/L43784_D20120125T2111/L43784_D20120125T211154.887Z_CS004_R000_tbb.h5'
output_file='/data/projects/COM_ALERT/pipeline/data/L43784_SB/L43784_D20120125T211154.887Z_CS004_sbtbb.h5'
subband_example_file='/data/projects/COM_ALERT/tbb/L597863_RS310_D20190809T090024.540Z_tbb.h5'
f_in=h5py.File(input_file,'r')
f_out=h5py.File(output_file,'w')
f_sample=h5py.File(subband_example_file,'r')

DEFAULT_ATTRS=[(u'TIME_RESOLUTION', 5.1200000000000001e-06),
 (u'GROUPTYPE', 'SubbandDataset'),
 (u'BANDWIDTH_UNIT', 'Hz'),
 (u'CENTRAL_FREQUENCY_UNIT', 'Hz'),
 (u'SAMPLES_PER_FRAME', 480),
 (u'TIME_RESOLUTION_UNIT', 's'),
 (u'BANDWIDTH', 195312.5)]

DIPOLE_ATTRS=[#u'ANTENNA_NORMAL_VECTOR',
 u'ANTENNA_POSITION',
# u'ANTENNA_POSITION_FRAME',
# u'ANTENNA_POSITION_UNIT',
 u'ANTENNA_ROTATION_MATRIX',
 u'GROUPTYPE',
 u'NOF_SUBBANDS',
 u'NYQUIST_ZONE',
 u'RCU_ID',
 u'RSP_ID',
 u'SAMPLE_FREQUENCY',
# u'SAMPLE_FREQUENCY_UNIT',
 u'STATION_ID']
# u'TILE_BEAM',
# u'TILE_BEAM_FRAME',
# u'TILE_BEAM_UNIT']

OBS_ATTRS=[u'DOC_VERSION',
 u'PROJECT_CO_I',
 u'FILETYPE',
 u'NOTES',
 u'DOC_NAME',
 u'FILENAME',
 u'OBSERVATION_NOF_STATIONS',
 u'OBSERVATION_FREQUENCY_MAX',
 u'OBSERVATION_END_MJD',
 u'FILTER_SELECTION',
 u'OBSERVATION_STATIONS_LIST',
 u'OBSERVATION_START_UTC',
 u'FILEDATE',
 u'TELESCOPE',
 u'ANTENNA_SET',
 u'OBSERVATION_END_UTC',
 u'PROJECT_PI',
 u'TARGETS',
 u'OBSERVATION_ID',
 u'OPERATING_MODE',
 u'PROJECT_CONTACT',
 u'SYSTEM_VERSION',
 u'OBSERVATION_NOF_BITS_PER_SAMPLE',
 u'GROUPTYPE',
 u'CLOCK_FREQUENCY',
 u'OBSERVATION_FREQUENCY_CENTER',
 u'PROJECT_ID',
 u'OBSERVATION_FREQUENCY_MIN',
 u'OBSERVATION_FREQUENCY_UNIT',
 u'CLOCK_FREQUENCY_UNIT',
 u'PROJECT_TITLE',
 u'OBSERVATION_START_MJD',
 u'NOF_STATIONS']

STAT_ATTRS=[u'GROUPTYPE',
 u'STATION_NAME',
 u'BEAM_DIRECTION',
 u'BEAM_DIRECTION_UNIT',
 u'BEAM_DIRECTION_FRAME',
 u'NOF_DIPOLES']


# FFT input
# Using a hanning window. For LOFAR a PPF would be more accurate.
window=np.hanning(1024)

obs_attrs=f_in.attrs
for key in OBS_ATTRS:
 if key in obs_attrs.keys():
      f_out.attrs.create(key,obs_attrs[key][0])
 else:
      print 'cannot obtain OBS',key
for k in f_in.keys():
    if 'Station' in k:
         #TODO: distinguish CS and RS
         K=k.upper().replace('N','N_CS')
         f_out.create_group(K)
         stat_attrs=f_in[k].attrs
         for key in STAT_ATTRS:
             if key in stat_attrs.keys():
                 f_out[K].attrs.create(key,stat_attrs[key][0])
             else:
                 print 'cannot obtain OBS',key

         for d in f_in[k].keys(): #TODO, run over all subbands
             D='DIPOLE_'+d
             f_out[K].create_group(D)
             timeseries_data=np.array(f_in[k][d])
             s=timeseries_data.shape
             dipole_attrs=f_in[k][d].attrs
             fftdata=np.fft.fft(timeseries_data.reshape(s[0]/1024,1024)*window)[:,512:]
             for key in DIPOLE_ATTRS:
                 if key in dipole_attrs.keys():
                      f_out[K][D].attrs.create(key,dipole_attrs[key])
                 elif key+'_VALUE' in dipole_attrs.keys():
                      f_out[K][D].attrs.create(key,dipole_attrs[key+'_VALUE'])
                 else:
                      print 'cannot obtain',key
             for SB in range(50,450): # TODO, expand range
                 S='SB_'+str(SB)
             	 f_out[K][D].create_dataset(S,data=fftdata[:,SB])     
                 for key,val in DEFAULT_ATTRS:
                        f_out[K][D][S].attrs.create(key,val)
                 # CENTRAL_FREQUENCY
                 key='CENTRAL_FREQUENCY'
                 val=1e8+SB*1e8/512
                 f_out[K][D][S].attrs.create(key,val)
                 # BAND NUMBER
                 key='BAND_NUMBER'
                 val=SB
                 # SLICE_NUMBER
                 key='SLICE_NUMBER'
                 val=dipole_attrs['SAMPLE_NUMBER']/1024
                 f_out[K][D][S].attrs.create(key,val) 
                 # TIME
                 key='TIME'
                 val=dipole_attrs['TIME']
                 f_out[K][D][S].attrs.create(key,val) 
                 
             


f_out.close()

