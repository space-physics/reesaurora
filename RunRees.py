#!/usr/bin/env python3
from dateutil.parser import parse
from numpy import array
from matplotlib.pyplot import show
import seaborn as sns
sns.color_palette(sns.color_palette("cubehelix"))
sns.set(context='talk', style='whitegrid',
        rc={'image.cmap': 'cubehelix_r'}) #for contour
#
from reesaurora.rees_model import reesiono,loadaltenergrid,energy_deg
from reesaurora.plots import fig11,plotA
from gridaurora.writeeigen import writeeigen
from gridaurora.solarangle import solarzenithangle
#
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

def runrees(t,glat,glon,isotropic,outfn,minalt,nalt,vlim):
    """
    inputs:
    t: time(s) to model, datetime() or ut1_unix


    """
    z,E = loadaltenergrid(minalt,nalt)

    q = reesiono(t, z, E, glat, glon,isotropic)

#%% outputs
    writeeigen(outfn,E,t,z,prates=q,tezs=None,latlon=(glat,glon))

    plotA(q,'Volume Production Rate {}  {} {}'.format(t,glat,glon),vlim)


def makefig11():
    from msise00.runmsis import rungtd1d

    z,E = loadaltenergrid(30.,200)

    E = array([1000,5000.])

    dens,temp = rungtd1d(t,altkm=z,glat=65.,glon=148.,
                         f107a=100.,f107=100.,ap=4.,
                     mass=48.,
                     tselecopts=array([1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],float)) #leave mass=48. !

#%% make figure 11 from Sergienko and Ivanov 1993
    Am,Lambda = energy_deg(E,isotropic,dens)

    fig11(E,Lambda)

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description='rees model of excitation rates for the aurora')
    p.add_argument('-t','--simtime',help='yyyy-mm-ddTHH:MM:SSZ time of sim',default='2013-01-01T12Z')
    p.add_argument('-c','--latlon',help='geodetic latitude/longitude (deg)',type=float,nargs=2,default=(65,-148))
    p.add_argument('--minalt',help='minimum altitude in grid [km]',type=float,default=30)
    p.add_argument('--nalt',help='Number of points in altitude grid',type=int,default=286)
    p.add_argument('-o','--outfn',help='give hdf5 filename to save eigenprofile production')
    p.add_argument('--isotropic',help='isotropic or non-isotropic pitch angle',action='store_true')
    p.add_argument('--vlim',help='plotting limits on energy dep and production plots',nargs=2,type=float,default=(1e-7,1e1))
    p = p.parse_args()

    t = parse(p.simtime)
#%%
    runrees(t,p.latlon[0],p.latlon[1], p.isotropic,
            p.outfn,p.minalt,p.nalt,p.vlim)

    print('solar zenith angle  {:.1f} '.format(solarzenithangle(p.simtime,p.latlon[0],p.latlon[1],0.)[0]))

    makefig11()

    show()
