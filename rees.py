#!/usr/bin/env python3
from __future__ import division
from os.path import expanduser
from numpy import logspace, loadtxt,tile
from dateutil.parser import parse
import h5py
#
from gridaurora.ztanh import setupz
from msise00.runmsis import rungtd1d
from rees_model import ionization_profile_from_flux

isotropic=False
"""
Discussion:
We need to know how incident flux from the magnetosphere ionizes various species
vs. altitude and inital beam energy in the ionosphere. The electron production rate q
may be related to auroral optical intensity via known proportions.

Variables defined:
q: e- production rate [cm^-3 s^-1]
Phi: "flux" [cm^-2 s^-1 eV^-1]
A: energy deposition matrix [eV cm^-1]

References:
Rees 1989
Wedlund et al "Electron Energy Spectra and Auroral Arcs" JGR 2013
Sergienko and Ivanov 1993 "A new approach to calculate the excitation of atmospheric gases by auroral electron impact"
"""
def reesmodel(dtime,altkm,E,glat,glon,f107a,f107,ap,mass,matfn,h5fn):

    dens,temp = rungtd1d(dtime,altkm,glat,glon,f107a,f107,ap,mass)
    rho = dens['Total'].values
    #%% Wedlund 2013 Eqn 4
    """
    each ij-th element of matrix A is defined in Eqn. 4
    """
    for i in range(altkm.size):
        for j,e in enumerate(E):
            A[i,j] = ( Lambda(s[i]/R[j]) * rho[i] * e *(1-albedo[j]) * dE[j] ) / ( R[e] * W[e])
    #%% Wedlund 2013 Eqn 1-3 (discretized form)
    q = A.dot(Phi)

def rees_aida(dtime,altkm,E,glat,glon,f107a,f107,ap,mass):
    dens,temp = rungtd1d(dtime,altkm,glat,glon,f107a,f107,ap,mass)
#%% temporary for testing with octave
#    from scipy.io import savemat
#    savemat('debug_msis.mat',{'z':altkm,'E':E,'O':dens['O'].values,
#                              'N2':dens['N2'].values,'O2':dens['O2'].values,
#                              'Total':dens['Total'].values},
#            oned_as='column')
#%% python port
    return ionization_profile_from_flux(E,dens,isotropic)

def plotA(A,z,E,ttxt):
    fg = figure()
    ax = fg.gca()
    hi = ax.imshow(A,extent=(E[0],E[-1],z[0],z[-1]),aspect='auto',origin='lower',
                                        interpolation='none')#,norm=LogNorm())
    fg.colorbar(hi,ax=ax)
    ax.set_title(ttxt)

    ax = figure().gca()
    ax.plot(A,z)
    ax.set_title(ttxt)

def writeA(A,z,E,dtime,ofn):
    with h5py.File(expanduser(ofn),'w',libver='latest') as f:
        f['/A']=A
        f['/z']=z
        f['/E'] =E
        f['/time']=str(dtime)
        f['/glat']=glat
        f['/glon']=glon

if __name__ == '__main__':
    from matplotlib.pyplot import figure,show
    #from matplotlib.colors import LogNorm

    import sys; sys.path.append('../transcarutils') #FIXME this approach is to avoid deeply nested transcarutils dependencies
    from readionoinit import getaltgrid #transcarutils
    import seaborn as sns
    sns.color_palette(sns.color_palette("cubehelix"))
    sns.set(context='talk', style='whitegrid',
            rc={'image.cmap': 'cubehelix_r'}) #for contour
#
    from argparse import ArgumentParser
    p = ArgumentParser(description='rees model of excitation rates for the aurora')
    p.add_argument('simtime',help='yyyy-mm-ddTHH:MM:SSZ time of sim',type=str,nargs='?',default='2013-03-31T09:00:00Z')
    p.add_argument('-c','--latlon',help='geodetic latitude/longitude (deg)',type=float,nargs=2,default=(65.12,-147.43 ))
    p.add_argument('--f107a',help=' 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)',type=float,default=107.6)
    p.add_argument('--f107',help='DAILY F10.7 FLUX FOR PREVIOUS DAY',type=float,default=126.0)
    p.add_argument('--ap',help='daily ap, 0-3hr, 3-6hr, 6-9hr, 9-12hr,12-33hr, 36-57hr',type=float,nargs=7,default=tile(63.7,7))
    p.add_argument('--mass',help=('MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS ' +
                       'CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL. '+
                         'MASS 17 IS Anomalous O ONLY.'),type=float,default=48)
    p.add_argument('--loadA',help='filename from already computed Rees A',type=str,default=None)
    p.add_argument('--loadE',help='loads E from transcar input',type=str,default='~/code/transcar/dir.transcar.server/BT_E1E2prev.csv')
    p.add_argument('--loadZ',help='loads Z from transcar input',type=str,default='~/code/transcar/dir.transcar.server/dir.input/conttanh.dat')
    p.add_argument('--h5',help='saves deposition matrix A as HDF5',type=str,default=None)
    p = p.parse_args()

    dtime = parse(p.simtime)
#%%
    if p.loadZ is not None:
        altkm = getaltgrid(expanduser(p.loadZ))
    else:
        altkm = setupz(286,90,1.5,11.1475)
    altkm = altkm[altkm<=1000] #keeps original spacing, but only auroral altitudes
#%% energy of beams
    if p.loadE is not None:
        E = loadtxt(expanduser(p.loadE),usecols=[0],delimiter=',')
    else:
        E = logspace(1.72,4.25,num=33,base=10)
#%%
    glat,glon = p.latlon
#%% from AIDA B. Gustavsson
    if p.loadA is not None:
        with h5py.File(expanduser(p.loadA),'r',libver='latest') as f:
            Aaida = f['/A'].value; zaida=f['/z'].value; Eaida=f['/E'].value
        plotA(Aaida, zaida, Eaida,'Rees deposition matrix')
    else: #compute A
        Aaida = rees_aida(dtime, altkm, E, glat, glon, p.f107a, p.f107, p.ap, p.mass)
        if p.h5:
            writeA(Aaida,altkm,E,dtime,p.h5)
        plotA(Aaida,altkm,E,'Rees deposition matrix')
#%% load JGR2013 data
    with h5py.File('../transcarutils/precompute/01Mar2011_FA.h5','r',libver='latest') as f:
        Ajgr = f['/Mp'].value; zjgr = f['/z'].value; Ejgr=f['/E'].value
    plotA(Ajgr.sum(axis=0),zjgr,Ejgr,'JGR2013 deposition matrix')
#%%
    #reesmodel(dtime, altkm, E, glat, glon, p.f107a, p.f107, p.ap, p.mass)

    show()
