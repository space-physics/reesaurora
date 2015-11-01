#!/usr/bin/env python3
from __future__ import division,absolute_import
from dateutil.parser import parse
from matplotlib.pyplot import show
import seaborn as sns
sns.color_palette(sns.color_palette("cubehelix"))
sns.set(context='talk', style='whitegrid',
        rc={'image.cmap': 'cubehelix_r'}) #for contour
#
from reesaurora.rees_model import reesiono,loadaltenergrid,plotA
from gridaurora.writeeigen import writeeigen
from gridaurora.solarangle import solarzenithangle
from gridaurora.readApF107 import readmonthlyApF107

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
    z,E = loadaltenergrid(minalt,nalt)

    q = reesiono(t, z, E, glat, glon, f107a, f107, ap, mass,isotropic)

#%% outputs
    writeeigen(outfn,E,t,z,prates=q,tezs=None,latlon=(glat,glon))

    plotA(q,'Rees deposition matrix',vlim)



if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description='rees model of excitation rates for the aurora')
    p.add_argument('simtime',help='yyyy-mm-ddTHH:MM:SSZ time of sim')
    p.add_argument('-c','--latlon',help='geodetic latitude/longitude (deg)',type=float,nargs=2,default=(65.12,-147.43 ))
#    p.add_argument('--f107a',help=' 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)',type=float)
#    p.add_argument('--f107',help='DAILY F10.7 FLUX FOR PREVIOUS DAY',type=float)
#    p.add_argument('--ap',help='daily ap, 0-3hr, 3-6hr, 6-9hr, 9-12hr,12-33hr, 36-57hr',type=float,nargs=7)
    p.add_argument('--mass',help=('MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS ' +
                       'CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL. '+
                         'MASS 17 IS Anomalous O ONLY.'),type=float,default=48)
    p.add_argument('--minalt',help='minimum altitude in grid [km]',type=float,default=30)
    p.add_argument('--nalt',help='Number of points in altitude grid',type=int,default=286)
    p.add_argument('-o','--outfn',help='give hdf5 filename to save eigenprofile production')
    p.add_argument('--isotropic',help='isotropic or non-isotropic pitch angle',action='store_true')
    p.add_argument('--vlim',help='plotting limits on energy dep and productio plots',nargs=2,type=float,default=(1e-7,1e1))
    p = p.parse_args()

    t = parse(p.simtime)
#%%
    f107Ap=readmonthlyApF107(int(str(t.year) + '{:02d}'.format(t.month)),'data/RecentIndices.txt')

    f107a = f107Ap['f107s']
    f107  = f107Ap['f107o']
    ap    = (f107Ap['Apo'],)*7

    runrees(t,p.latlon[0],p.latlon[1],
            f107a,f107,ap,p.mass,
            p.isotropic,
            p.outfn,p.minalt,p.nalt,p.vlim)

    print('solar zenith angle  {:.1f}  f107a {}  f107 {}  Ap {}'.format(solarzenithangle(
           p.simtime,p.latlon[0],p.latlon[1],0.)[0],
            f107a,f107,ap))

    show()
