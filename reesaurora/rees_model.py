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
from pandas import Panel
from numpy import (gradient,array,linspace,zeros,diff,append,empty,arange,log10,exp,nan,
                   logspace,atleast_1d,ndarray,copy)
from scipy.interpolate import interp1d
#
from gridaurora.ztanh import setupz
from msise00.runmsis import rungtd1d
from gridaurora.readApF107 import readmonthlyApF107
try:
    from glowaurora.runglow import glowalt
except ImportError as e:
    logging.error(e)

def reesiono(T,altkm:ndarray,E:ndarray,glat:float,glon:float,isotropic:bool):
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

    qPanel = Panel(items=T,
                   major_axis=E,
                   minor_axis=altkm)
#%% loop
    for t in T:
        f107Ap=readmonthlyApF107(t)
        f107a = f107Ap['f107s']
        f107  = f107Ap['f107o']
        ap    = (f107Ap['Apo'],)*7

        dens,temp = rungtd1d(t,altkm,glat,glon,f107a,f107,ap,
                             mass=48.,
                             tselecopts=array([1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],float)) #leave mass=48. !

        q = ionization_profile_from_flux(E,dens,isotropic)
        qPanel.loc[t,:,:] = q.T

    return qPanel

#TODO check that "isotropic" is handled consistently with original code
def ionization_profile_from_flux(E,dens,isotropic):
    """
    simple model for volume emission as function of altitude.
    After Sergienko and Ivanov 1993 and Gustavsson AIDA_TOOLs
    """

    if ((E<1e2) | (E>1e4)).any():
        logging.warning('Sergienko & Ivanov 1993 covered E \in [100,10000] eV')

    if (dens.index>700.).any():
        logging.warning('Sergienko & Ivanov 1993 assumed electron source was at altitude 700km.')

#%% Table 1 Sergienko & Ivanov 1993, rightmost column
    # mean energy per ion-electron pair
    E_cost_ion = array([36.8,26.8,28.2]) # N_2, O, O_2

    ki = array([1, 0.7, 0.4])

    dE = diff(E); dE = append(dE,dE[-1])

    Partitioning = partition(dens,ki)

#%% Calculate the energy deposition as a function of altitude
    qE = empty((dens.shape[0],E.size)) # Nalt x Nenergy
    for i,(e,d) in enumerate(zip(E,dE)):
        Ebins = linspace(e,e+d,20) #make a subset of fine resolution energy bins within bigger  energy bins
        #for isotropic or field aligned electron beams
        Am = energy_deg(Ebins,isotropic,dens) # Nsubenergy x Naltitude

        q= Am.sum(axis=0) #sum over the interim energy sub-bins
        q *= (Partitioning/E_cost_ion).sum(axis=1) #effect of ion chemistry at each altitude
        qE[:,i] = q
    return qE

def energy_deg(E,isotropic,dens):
    """
    energy degradation of precipitating electrons
    """
    #atmp = DataFrame(index=dens.index)
    #atmp[['N2','O','O2']] = dens[['N2','O','O2']]/1e6
    atmp = dens['Total'].values/1e3

    N_alt0 = atmp.shape[0]
    zetm = zeros(N_alt0)
    dH = gradient(dens.index)
    for i in range(N_alt0-1,0,-1): #careful with these indices!
        dzetm = (atmp[i] +atmp[i-1])*dH[i-1]*1e5/2
        zetm[i-1] = zetm[i] + dzetm

    alb = albedo(E,isotropic)

    Am = zeros((E.size,N_alt0))
    D_en = gradient(E)
    r = PitchAngle_range(E,isotropic)

    chi = zetm / r[:,None]

    Lambda = lambda_comp(chi,E,isotropic)[0]

    Am = atmp * Lambda * E[:,None] * (1-alb[:,None])/r[:,None]

    Am[0,:] *= D_en[0]/2.
    Am[-1,:]*= D_en[-1]/2.
    Am[1:-2,:] *= (D_en[1:-2]+D_en[0:-3])[:,None]/2.
    return Am

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

def albedo(E,isotropic,fn='data/SergienkoIvanov.h5'):

    with h5py.File(fn,'r',libver='latest') as h:
        albedoflux = h['/albedo/flux'][int(isotropic),:]
        LE = h['/albedo/E']

        Emax = 10**LE[-1]
        E = copy(E)
        E[E>Emax] = Emax

        falb=interp1d(LE, albedoflux, kind='linear',bounds_error=False,fill_value=nan)
        alb = falb(log10(E))

    return alb

def lambda_comp(chi,E,isotropic,fn='data/SergienkoIvanov.h5'):
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

def partition(dens,ki):
    P = ki[[0,2,1]]*dens[['N2','O','O2']].values
    return P / P.sum(axis=1)[:,None]

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

    z = z[z <= 1000] #keeps original spacing, but only auroral altitudes
#%% energy of beams
    if special_grid.lower()=='transcar':
        E = logspace(1.72,4.25,num=33,base=10)
    else:
        E = logspace(1.72,6.,num=81,base=10)

    return z,E