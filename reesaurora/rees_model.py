#!/usr/bin/env python3
"""
 ionization_profiles_from_flux - simple model for volume emission as function of altitude.
   After Sergienko and Ivanov 1993
   a massively speeded up implementation after the AIDA_TOOLS package by Gustavsson, Brandstrom, et al
"""
import logging
import h5py
from dateutil.parser import parse
from datetime import datetime
from xarray import DataArray
from numpy import (gradient,array,linspace,zeros,diff,empty,append,log10,exp,nan,
                   logspace,atleast_1d,ndarray,copy)
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
#
from gridaurora.ztanh import setupz
from msise00.runmsis import rungtd1d
from gridaurora.readApF107 import readmonthlyApF107
try:
    from glowaurora.runglow import glowalt
except ImportError as e:
    logging.error(e)

from .plots import fig11

species =['N2','O','O2']

def reesiono(T,altkm:ndarray,E:ndarray,glat:float,glon:float,isotropic:bool,verbose:int,datfn):
    #other assertions covered inside modules
    assert isinstance(isotropic,bool)

    if isinstance(T,str):
        T=parse(T)
    T = atleast_1d(T)
    assert isinstance(T[0],datetime)
#%% MSIS
    if isotropic:
        logging.debug('isotropic pitch angle flux')
    else:
        logging.debug('field-aligned pitch angle flux')

    Qt = DataArray(data=empty((T.size,len(species),altkm.size,E.size)),
                   coords=[T,species,altkm,E],
                   dims=['time','species','altkm','energy'])
#%% loop
    for t in T:
        f107Ap=readmonthlyApF107(t)
        f107a = f107Ap['f107s']
        f107  = f107Ap['f107o']
        ap    = (f107Ap['Apo'],)*7

        dens,temp = rungtd1d(t,altkm,glat,glon,f107a,f107,ap,
                             mass=48.,
                             tselecopts=array([1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],float)) #leave mass=48. !

        Q = ionization_profile_from_flux(E,dens,isotropic,species,datfn,verbose)
        Qt.loc[t,...] = Q

    return Qt

def ionization_profile_from_flux(E,dens,isotropic,species,datfn,verbose):
    """
    simple model for volume emission as function of altitude.
    After Sergienko and Ivanov 1993 and Gustavsson AIDA_TOOLs
    """
    if ((E<100) | (E>1e4)).any():
        logging.warning('Sergienko & Ivanov 1993 covered E \in [100,10000] eV')

    if (dens.altkm>700.).any():
        logging.warning('Sergienko & Ivanov 1993 assumed electron source was at altitude 700km.')

#%% Table 1 Sergienko & Ivanov 1993, rightmost column
    # mean energy per ion-electron pair
    E_cost_ion = {'N2':36.8,'O2':28.2,'O':26.8}
#%% Eqn 7, Figure 6
    k = {'N2':1.,'O2':0.7,'O':0.4}

    dE = diff(E); dE = append(dE,dE[-1])

    Peps = partition(dens,k,E_cost_ion) #Nalt x Nspecies

    # NOTE kg m^-2 *10 = g cm^-2
    rho = dens.loc[:,'Total'] * 10. # mass density of atmosphere [g cm^-3]
#%% Calculate the energy deposition "loss" as a function of altitude
    """
    Implement through Eqn 8.
    Q: per species per state
    """
    Q = DataArray(data = empty((len(species),dens.shape[0],E.size)),
                  coords=[species, dens.altkm, E],
                  dims=['species','altkm','energy']) # Nalt x Nenergy x Nspecies

    for i,(e,d) in enumerate(zip(E,dE)):
        Ebins = linspace(e,e+d,20) #make a subset of fine resolution energy bins within bigger  energy bins

        #for isotropic or field aligned electron beams
        # Eqn 6,A1
        W = energy_deg(Ebins,isotropic,rho,datfn,verbose) # Nsubenergy x Naltitude

        #Eqn A4
        for s in species:
            # we sum W across the subenergies (numerical integration to help large energy grid step size)
            Q.loc[s,:,e] = Peps.loc[:,s] * rho * W.sum(axis=0) # production rate [cm^-3 s^-1] due to gas mass density impacts at each altitude

    return Q

def energy_deg(E,isotropic,rho,datfn,verbose):
    """
    energy degradation of precipitating electrons -- kernel of Eqn A4,6 implementing A1
    rho: MASS density g cm^-3
    """

    z = rho.altkm.values
    #rho=rho.values #so Numpy can cope

    dE = gradient(E)

    Rng = PitchAngle_range(E,isotropic) # Eqn A3, Table 7
    alb = albedo(E,isotropic,datfn)

    #scattering depth (Rees  Physics & Chem of Upper Atm. Eqn 3.3.5 pg. 40)
    zetm = -cumtrapz(rho[::-1],z[::-1],initial=-0.)[::-1]
#%% Eqn 8 as Eqn A4
    chi = zetm / Rng[:,None]
    Lambda = lambda_comp(chi,E,isotropic,datfn)[0]
    if verbose > 1: #dozens of plots (Nenergy)
        if isotropic:
            Lambda_i=[Lambda[0,:]]; Lambda_m=None
        else:
            Lambda_m=[Lambda[0,:]]; Lambda_i=None
        fig11([E[0]],[chi[0,:]],Lambda_m,Lambda_i,None)
#%%  Eqn A1
    #Nsubenergy x Nalt
    W = Lambda * E[:,None] * (1-alb[:,None]) / Rng[:,None]

    W[0,:] *= dE[0]/2. #lowest subenergy
    W[-1,:]*= dE[-1]/2.
    W[1:-2,:] *= (dE[1:-2] + dE[0:-3])[:,None] / 2.

    return W

def PitchAngle_range(E,isotropic):
    """
    Eqn A3 & Table 7 from Sergienko & Ivanov 1993
    Note E is keV in Eqn A3, hence the divide by 1000

    output:
    range(E) [g/cm^2]
    """
    kE = E/1000.

    B2 = 9.48e-2
    B3 = -1.57

    # paper says this, but doesn't match Fig. 13. Bjorn had what's used.
    B1=1.804e-6 if isotropic else 2.16e-6

    #B1 = 1.64e-6 if isotropic else 2.16e-6

    return B1*kE**1.67 * (1. + B2*kE**B3)

    #pr= 1.64e-6 if isotropic else 2.16e-6
    #return pr * (E/1e3)**1.67 * (1 + 9.48e-2 * E**-1.57)

def albedo(E,isotropic,fn):

    with h5py.File(str(fn),'r',libver='latest') as h:
        albedoflux = h['/albedo/flux'][int(isotropic),:]
        LE = h['/albedo/E']

        Emax = 10**LE[-1]
        E = copy(E)
        E[E>Emax] = Emax

        falb=interp1d(LE, albedoflux, kind='linear',bounds_error=False,fill_value=nan)
        alb = falb(log10(E))

    return alb

def lambda_comp(chi,E,isotropic,fn):
    """
    Implements Eqn. A2 from Sergienko & Ivanov 1993

    field-aligned monodirectional: interpolated over energies from 48.9 eV to 5012 eV

    Isotropic: interpolated over 48.9ev to 1000 eV

    "Param_m" and "Param_i" are from Table 6 of Sergienko & Ivanov 1993

    """
    with h5py.File(str(fn),'r',libver='latest') as h:
#%% choose isotropic or monodirectional
        if isotropic:
            P = h['isotropic/C']
            LE =h['isotropic/E']
        else:
            P = h['monodirectional/C']
            LE =h['monodirectional/E']
#%% more robust way to handle too-high values, like paper appears to do
        Emax = 10**LE[-1]
        E = copy(E)
        E[E>Emax] = Emax
#%% interpolate  -- use NaN as a sentinal value
        fC=interp1d(LE,P,kind='linear',axis=1,bounds_error=False,fill_value=nan)
        C = fC(log10(E))
        """
        the section below finally implements Eqn. A2 from the Sergienko & Ivanov 1993 paper.
        We create a plot mimicing Fig. 11 from this paper.
        """
        lam = ((C[0,:][:,None]*chi + C[1,:][:,None]) *
                exp(C[2,:][:,None]*chi**2 + C[3,:][:,None]*chi))

    return lam,C

def partition(dens,k,cost):
    """
    Implement Eqn 7 Sergienko 1993

    N: NUMBER density [cm^-3] vs. altitude for each species
    k: correction factors vs. Monte Carlo for Sergienko 1993
    cost: energization cost

    output:
    P_i(h) / epsilon
    """
    # m^-3 /1e6 = cm^-3
    N = dens.loc[:,species]/1e6  #[cm^-3]

    num=DataArray(data=empty((N.shape[0],len(species))),
                  coords=[N.altkm,species],
                  dims=['altkm','species'])
    for i in species:
        num.loc[:,i] = k[i]*N.loc[:,i]

    den=num.sum('species')

    P = num.copy()
    for i in species:
        P.loc[:,i] = num.loc[:,i]/den
#%% energization cost
    Peps = num.copy()
    for i in species:
        Peps.loc[:,i] = P.loc[:,i] / cost[i]

    return Peps

def loadaltenergrid(minalt=90,Nalt=286,special_grid=''):
    """
    makes a tanh-spaced grid (see setupz for info)

    minalt: [km] minimum altiude in grid (e.g. 90)
    Nalt: number of points in grid
    special_grid: use same grid as 'transcar' or 'glow'
    """
    assert isinstance(special_grid,(str,None))
    #%% altitude
    if special_grid.lower()=='transcar':
        z = setupz(286,90,1.5,11.1475)
    elif special_grid.lower()=='glow':
        z = glowalt()
    else:
        z = setupz(Nalt,minalt,1.5,11.1475)

    z = z[z <= 700] #keeps original spacing, but only auroral altitudes to source at z=700km
#%% energy of beams
    if special_grid.lower()=='transcar':
        E = logspace(1.72,4.25,num=33,base=10)
    else:
        E = logspace(1.69,4.,num=81,base=10)

    return z,E