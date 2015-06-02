#!/usr/bin/env python3
"""
 ionization_profiles_from_flux - simple model for volume emission as function of altitude.
   After Sergienko and Ivanov 1993
   a massively speeded up implementation after the AIDA_TOOLS package by Gustavsson, Brandstrom, et al
"""
from pandas import DataFrame
from numpy import gradient,array,linspace,zeros,diff,append,empty,arange,log10,exp,nan
from scipy.interpolate import interp1d
#
from oct2py import Oct2Py
oc = Oct2Py(executable='octave',oned_as='column',convert_to_float=True,
                                        timeout=30)


def ionization_profile_from_flux(E,dens,isotropic):
    E_cost_ion = array([36.8,26.8,28.2])
    ki = array([1, 0.7, 0.4])

    dE = diff(E); dE = append(dE,dE[-1])

    Partitioning = partition(dens,ki)

#%% First calculate the energy deposition as a function of altitude
    qall = empty((dens.shape[0],E.size))
    for i,(e,d) in enumerate(zip(E,dE)):
        Ebins = linspace(e,e+d,20)
        #for both isotropic and field aligned electron beams
        Am = energy_deg(Ebins,isotropic,dens)

        q= Am.sum(axis=0)
        q *= (Partitioning/E_cost_ion).sum(axis=1)
        qall[:,i] = q
    return qall

def energy_deg(E,isotropic,dens):
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
    r = oc.Pat_range(E,isotropic+1)

    hi = zetm/r

    Lambda = lambda_comp(hi,E,isotropic)

    Am = atmp.iloc[:,-1].values * Lambda * E[:,None] * (1-alb[:,None])/r

#    for i in range(N_alt0):
#        Am[:,i] = atmp.iat[i,-1] * Lambda[:,i] * E[:,None] * (1-alb)/r

    Am[0,:] *= D_en[0]/2
    Am[-1,:]*= D_en[-1]/2
    Am[1:-2,:] *= (D_en[1:-2]+D_en[0:-3])[:,None]/2
    return Am

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
    logE_m= append(1.69,arange(1.8,3.7+0.1,0.1))
    Param_m =array(
        [[1.43,1.51,1.58,1.62,1.51,1.54,1.18,1.02, 0.85, 0.69, 0.52,0.35,0.21,0.104,0.065,0.05,0.04,0.03,0.03, 0.025,0.021],
         [0.83,0.77,0.72,0.67,0.63,0.59,0.56,0.525,0.495,0.465,0.44,0.42,0.40,0.386,0.37, 0.36,0.35,0.34,0.335,0.325,0.32],
         [-0.025,-0.030, -0.040, -0.067, -0.105, -0.155, -0.210, -0.275, -0.36, -0.445, -0.51,-0.61, -0.69, -0.77, -0.83, -0.865, -0.90,-0.92, -0.935, -0.958, -0.96],
         [-1.67,-1.65,-1.62,-1.56,-1.46,-1.35,-1.20,-0.98,-0.70,-0.37,-0.063,0.39,0.62,0.92,1.11,1.25,1.36,1.44,1.50,1.55,1.56]])

    logE_i=append(1.69, arange(1.8,3.0+0.1,0.1))
    Param_i =array(
        [[0.041, 0.051, 0.0615, 0.071, 0.081, 0.09, 0.099, 0.1075, 0.116, 0.113, 0.13, 0.136, 0.139, 0.142],
         [1.07, 1.01, 0.965, 0.9, 0.845, 0.805, 0.77, 0.735, 0.71, 0.69, 0.67, 0.665, 0.66, 0.657],
         [-0.064, -0.1, -0.132, -0.171, -0.2, -0.221, -0.238, -0.252, -0.261, -0.267, -0.271, -0.274, -0.276, -0.277],
         [-1.054, -0.95, -0.845, -0.72, -0.63, -0.54, -0.475, -0.425, -0.38, -0.345, -0.319, -0.295, -0.28, -0.268]])

    logE=log10(E)

    if isotropic:
        P = Param_i
        LE = logE_i
        Emax = 1000.
    else:
        P = Param_m
        LE = logE_m
        Emax = 5000.

    fC=interp1d(LE,P,kind='linear',axis=1,bounds_error=False,fill_value=nan)
    C = fC(logE)

    lam = ((C[0,:][:,None]*hi + C[1,:][:,None]) *
            exp(C[2,:][:,None]*hi**2 + C[3,:][:,None]*hi))

    badind = E>Emax
    lam[badind] = ((P[0,-1]*hi[badind] + P[1,-1]) *
            exp(P[2,-1]*hi[badind]**2 +P[3,-1]*hi[badind]))
    return lam

def partition(dens,ki):
    k = array([ki[0],ki[2],ki[1]])
    P = k*dens[['N2','O','O2']].values
    return P / P.sum(axis=1)[:,None]