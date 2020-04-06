# TBBs
LOFAR TBB commissioning

## TBB beamforming individual station
Beamforming TBB data for station CS001, polarization 0, HBA0:
```
python bf_data.py /data/projects/COM_ALERT/tbb/L597863_CS001_D20200317T143637.202Z*tbb.h5 -p 0 -s HBA0
```
It will produce a tbb beamformed .h5 file and the dynamic spectrum in .beam (pycrtools) format.
