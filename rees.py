#!/usr/bin/env python3
from __future__ import division,absolute_import
from numpy import tile
from matplotlib.pyplot import show
import seaborn as sns
sns.color_palette(sns.color_palette("cubehelix"))
sns.set(context='talk', style='whitegrid',
        rc={'image.cmap': 'cubehelix_r'}) #for contour
#
from rees_model import reesiono,loadaltenergrid,plotA
from gridaurora.writeeigen import writeeigen

isotropic=False
"""
Models unit incident flux from the magnetosphere ionizing N2,O,O2
vs. altitude and inital beam energy in the ionosphere.
The electron production rate q is related to auroral optical intensity via
known proportions.

Variables:
q: e- production rate [cm^-3 s^-1]
Phi: "flux" [cm^-2 s^-1 eV^-1]
A: energy deposition matrix [eV cm^-1]

References:
Rees 1989
Wedlund et al "Electron Energy Spectra and Auroral Arcs" JGR 2013
Sergienko and Ivanov 1993 "A new approach to calculate the excitation of atmospheric gases by auroral electron impact"
"""

def runrees(t,glat,glon,f107a,f107,ap,mass,isotropic,outfn,minalt,nalt,vlim):
    """
    inputs:
    t: time(s) to model, datetime() or ut1_unix


    """
    altkm,E = loadaltenergrid(minalt,nalt)

    Aaida = reesiono(t, altkm, E, glat, glon, f107a, f107, ap, mass,isotropic)

#%% outputs
    writeeigen(outfn,E,t,tezs=None,latlon=(glat,glon))

    plotA(Aaida,altkm,E,'Rees deposition matrix',vlim)



if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description='rees model of excitation rates for the aurora')
    p.add_argument('simtime',help='yyyy-mm-ddTHH:MM:SSZ time of sim',nargs='?',default='2013-03-31T09:00:00Z')
    p.add_argument('-c','--latlon',help='geodetic latitude/longitude (deg)',type=float,nargs=2,default=(65.12,-147.43 ))
    p.add_argument('--f107a',help=' 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)',type=float,default=107.6)
    p.add_argument('--f107',help='DAILY F10.7 FLUX FOR PREVIOUS DAY',type=float,default=126.0)
    p.add_argument('--ap',help='daily ap, 0-3hr, 3-6hr, 6-9hr, 9-12hr,12-33hr, 36-57hr',type=float,nargs=7,default=tile(63.7,7))
    p.add_argument('--mass',help=('MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS ' +
                       'CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL. '+
                         'MASS 17 IS Anomalous O ONLY.'),type=float,default=48)
    p.add_argument('--minalt',help='minimum altitude in grid [km]',type=float,default=30)
    p.add_argument('--nalt',help='Number of points in altitude grid',type=int,default=286)
    p.add_argument('-o','--outfn',help='give hdf5 filename to save eigenprofile production')
    p.add_argument('--isotropic',help='isotropic or non-isotropic pitch angle',action='store_true')
    p.add_argument('--vlim',help='plotting limits on energy dep and productio plots',nargs=2,type=float,default=(1e-7,1e1))
    p = p.parse_args()

    runrees(p.simtime,p.latlon[0],p.latlon[1],
            p.f107a,p.f107,p.ap,p.mass,p.isotropic,
            p.outfn,p.minalt,p.nalt,p.vlim)

    show()
