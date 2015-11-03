#!/usr/bin/env python3
"""
 ionization_profiles_from_flux - simple model for volume emission as function of altitude.
   After Sergienko and Ivanov 1993
   a massively speeded up implementation after the AIDA_TOOLS package by Gustavsson, Brandstrom, et al
"""
from __future__ import division,absolute_import
from six import string_types
from pandas import DataFrame,Panel
from numpy import (gradient,array,linspace,zeros,diff,append,empty,arange,log10,exp,nan,
                   logspace,atleast_1d)
from scipy.interpolate import interp1d
from matplotlib.pyplot import figure
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator
#
from gridaurora.ztanh import setupz
from msise00.runmsis import rungtd1d
from gridaurora.readApF107 import readmonthlyApF107
try:
    from glowaurora.runglow import glowalt
except:
    pass

def reesiono(T,altkm,E,glat,glon,isotropic):
    #other assertions covered inside modules
    assert isinstance(isotropic,bool)
    T = atleast_1d(T)
#%% MSIS
    if isotropic:
        print('isotropic pitch angle flux')
    else:
        print('field-aligned pitch angle flux')

    qPanel = Panel(items=T,
                   major_axis=E,
                   minor_axis=altkm)
#%% loop
    for t in T:
        f107Ap=readmonthlyApF107(t,'data/RecentIndices.txt')
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
    E_cost_ion = array([36.8,26.8,28.2])
    ki = array([1, 0.7, 0.4])

    dE = diff(E); dE = append(dE,dE[-1])

    Partitioning = partition(dens,ki)

#%% First calculate the energy deposition as a function of altitude
    qE = empty((dens.shape[0],E.size)) # Nalt x Nenergy
    for i,(e,d) in enumerate(zip(E,dE)):
        Ebins = linspace(e,e+d,20)
        #for isotropic or field aligned electron beams
        Am = energy_deg(Ebins,isotropic,dens)

        q= Am.sum(axis=0) #sum over the interim energy sub-bins
        q *= (Partitioning/E_cost_ion).sum(axis=1) #effect of ion chemistry at each altitude
        qE[:,i] = q
    return qE

def energy_deg(E,isotropic,dens):
    """
    energy degradation of precipitating electrons
    """
    atmp = DataFrame(index=dens.index)
    atmp[['N2','O','O2']] = dens[['N2','O','O2']]/1e6
    atmp['Total'] = dens['Total']/1e3

    N_alt0 = atmp.shape[0]
    zetm = zeros(N_alt0)
    dH = gradient(atmp.index)
    for i in range(N_alt0-1,0,-1): #careful with these indices!
        dzetm = (atmp.iat[i,-1]+atmp.iat[i-1,-1])*dH[i-1]*1e5/2
        zetm[i-1] = zetm[i] + dzetm

    alb = albedo(E,isotropic)

    Am = zeros((E.size,N_alt0))
    D_en = gradient(E)
    r = Pat_range(E,isotropic)

    hi = zetm / r[:,None]

    Lambda = lambda_comp(hi,E,isotropic)

    Am = atmp.iloc[:,-1].values * Lambda * E[:,None] * (1-alb[:,None])/r[:,None]

#    for i in range(N_alt0):
#        Am[:,i] = atmp.iat[i,-1] * Lambda[:,i] * E[:,None] * (1-alb)/r

    Am[0,:] *= D_en[0]/2.
    Am[-1,:]*= D_en[-1]/2.
    Am[1:-2,:] *= (D_en[1:-2]+D_en[0:-3])[:,None]/2.
    return Am

def Pat_range(E,isotropic):
    pr= 1.64e-6 if isotropic else 2.16e-6
    return pr * (E/1e3)**1.67 * (1 + 9.48e-2 * E**-1.57)

def albedo(E,isotropic):
    isotropic = int(isotropic)
    logE_p=append(1.69, arange(1.8,3.7+0.1,0.1))
    Param=array(
      [[0.352, 0.344, 0.334, 0.320, 0.300, 0.280, 0.260, 0.238, 0.218, 0.198, 0.180, 0.160, 0.143, 0.127, 0.119, 0.113, 0.108, 0.104, 0.102, 0.101, 0.100],
       [0.500, 0.492, 0.484, 0.473, 0.463, 0.453, 0.443, 0.433, 0.423, 0.413, 0.403, 0.395, 0.388, 0.379, 0.378, 0.377, 0.377, 0.377, 0.377, 0.377, 0.377]])
    logE=log10(E)

    falb=interp1d(logE_p,Param[isotropic,:],kind='linear',bounds_error=False,fill_value=nan)
    alb = falb(logE)
    alb[logE>logE_p[-1]] = Param[isotropic,-1]

    return alb

def lambda_comp(hi,E,isotropic):
    """
    interpolated over energies from 48.9 eV to 5012 eV
    for isotropic and field-aligned precipitation
    """
#%% field-aligned
    logE_m= append(1.69,arange(1.8,3.7+0.1,0.1))
    Param_m =array(
        [[1.43,1.51,1.58,1.62,1.51,1.54,1.18,1.02, 0.85, 0.69, 0.52,0.35,0.21,0.104,0.065,0.05,0.04,0.03,0.03, 0.025,0.021],
         [0.83,0.77,0.72,0.67,0.63,0.59,0.56,0.525,0.495,0.465,0.44,0.42,0.40,0.386,0.37, 0.36,0.35,0.34,0.335,0.325,0.32],
         [-0.025,-0.030, -0.040, -0.067, -0.105, -0.155, -0.210, -0.275, -0.36, -0.445, -0.51,-0.61, -0.69, -0.77, -0.83, -0.865, -0.90,-0.92, -0.935, -0.958, -0.96],
         [-1.67,-1.65,-1.62,-1.56,-1.46,-1.35,-1.20,-0.98,-0.70,-0.37,-0.063,0.39,0.62,0.92,1.11,1.25,1.36,1.44,1.50,1.55,1.56]]
         )
#%% isotropic
    """
        interpolated over energies from 48.9 eV to 1000 eV
    """
    logE_i=append(1.69, arange(1.8,3.0+0.1,0.1))
    Param_i =array(
        [[0.041, 0.051, 0.0615, 0.071, 0.081, 0.09, 0.099, 0.1075, 0.116, 0.113, 0.13, 0.136, 0.139, 0.142],
         [1.07, 1.01, 0.965, 0.9, 0.845, 0.805, 0.77, 0.735, 0.71, 0.69, 0.67, 0.665, 0.66, 0.657],
         [-0.064, -0.1, -0.132, -0.171, -0.2, -0.221, -0.238, -0.252, -0.261, -0.267, -0.271, -0.274, -0.276, -0.277],
         [-1.054, -0.95, -0.845, -0.72, -0.63, -0.54, -0.475, -0.425, -0.38, -0.345, -0.319, -0.295, -0.28, -0.268]]
         )

    logE=log10(E)

    if isotropic:
        P = Param_i
        LE = logE_i
        Emax = 1000.
    else:
        P = Param_m
        LE = logE_m
        Emax = 5000.
#%% interpolate
    fC=interp1d(LE,P,kind='linear',axis=1,bounds_error=False,fill_value=nan)
    C = fC(logE)
#%% low energy
    lam = ((C[0,:][:,None]*hi + C[1,:][:,None]) *
            exp(C[2,:][:,None]*hi**2 + C[3,:][:,None]*hi))
#%% high energy
    badind = E>Emax
    lam[badind] = (
                   (P[0,-1]*hi[badind] + P[1,-1]) *
                   exp(P[2,-1]*hi[badind]**2 + P[3,-1]*hi[badind])
                   )

    return lam

def partition(dens,ki):
    P = ki[[0,2,1]]*dens[['N2','O','O2']].values
    return P / P.sum(axis=1)[:,None]

def plotA(q,ttxt,vlim):
    E=q.major_axis.values
    z=q.minor_axis.values
    Q=q.values.squeeze().T
#%%
    def _doax(ax):
       ax.yaxis.set_major_locator(MultipleLocator(100))
       ax.yaxis.set_minor_locator(MultipleLocator(20))
       ax.set_xscale('log')
       ax.set_ylabel('altitude [km]')
       ax.set_title(ttxt)

    fg = figure()
    ax = fg.gca()
    hi = ax.pcolormesh(E,z,Q,
                       vmin=vlim[0],vmax=vlim[1],
                       norm=LogNorm())
    c=fg.colorbar(hi,ax=ax)
    c.set_label('Volume Production Rate')
    ax.set_xlabel('beam energy [eV]')
    ax.autoscale(True,tight=True) #fill axes
    _doax(ax)
#%% same data, differnt plot
    ax = figure().gca()
    ax.plot(Q,z)
    ax.set_xlim(vlim)
    ax.set_xlabel('Energy Deposition')
    _doax(ax)

def loadaltenergrid(minalt=90,Nalt=286,special_grid=''):
    """
    makes a tanh-spaced grid (see setupz for info)

    minalt: [km] minimum altiude in grid (e.g. 90)
    Nalt: number of points in grid
    special_grid: use same grid as 'transcar' or 'glow'
    """
    assert isinstance(special_grid,string_types)
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