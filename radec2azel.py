#! /usr/bin/env python
import pycrtools as cr
from pycrtools import tools
from pytmf import *
import numpy as np
from pycrtools import metadata as md

# From pipelines/frats_event
def getRADEC_2_AZEL(alpha, delta, utc, angl_offset=0.0):
    ''' Quick funciton to change coordinates from RADEC to AZEL.
    '''

    phi = deg2rad(52.915122495) #(LOFAR Superterp)
    L = deg2rad(6.869837540)

    angl_offset = deg2rad(angl_offset)
    alpha += angl_offset
    #delta = delta+angl_offset

    #----------------------------------------------------
    # Make input and output arrays for conversion
    equatorial = cr.hArray([alpha, delta])
    horizontal = equatorial.new()

    #----------------------------------------------------
    # Convert all coordinates in the input array
    # assumes difference between UTC and UT is 0 (hence ut1_utc=0.)
    cr.hEquatorial2Horizontal(horizontal, equatorial, utc, 0., L, phi)

    if horizontal[0] < 0:
        horizontal[0] += 2. * np.pi # Need to check the definitions used for positive azimut angles.

    azel = [horizontal[0],horizontal[1]]
    pointings = [{'az': azel[0], 'el': azel[1]}]

    return pointings

def azel2beamweights(azel,station,freq,antennaset,lofarcentered=False):
    #antennaset = HBA0, HBA1 or HBA_DUAL
    if lofarcentered:
        antpos=md.getLofarCenterRelativeAntennaPosition(station,antennaset,True)
    else:
        antpos=md.getStationCenterRelativeAntennaPosition(station,antennaset,True)
    delays=cr.hArray(float, [96])
    pointings=cr.hArray([cr.convert(coords, "CARTESIAN") for coords in azel][0])
    cr.hGeometricDelays(delays, antpos, pointings, True)
    phases=cr.hArray(float,[96])
    phases.delaytophase(cr.hArray([freq]),delays)
    weights=cr.hArray(complex, 96)
    weights.phasetocomplex(phases)
    return weights 

#t=tools.strdate2jd('2019-11-26 15:22:00')
#ra=hms2rad(3,32,59.34)
#dec=dms2rad(54,34,43.2)
#azel=getRADEC_2_AZEL(ra,dec,t)
#print ra,dec,t,azel
#station="RS106"
#substation='HBA'
#freq=185000000.0#400*0.1953125+100e6
#weights=azel2beamweights(azel,station,freq,substation)
