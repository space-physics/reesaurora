#!/usr/bin/env python

isotropic = False
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
# def reesmodel(dtime,altkm,E,glat,glon,f107a,f107,ap,mass,matfn,h5fn):
#
#    dens,temp = rungtd1d(dtime,altkm,glat,glon,f107a,f107,ap,mass)
#    rho = dens['Total'].values
#    #%% Wedlund 2013 Eqn 4
#    """
#    each ij-th element of matrix A is defined in Eqn. 4
#    """
#    for i in range(altkm.size):
#        for j,e in enumerate(E):
#            A[i,j] = ( Lambda(s[i]/R[j]) * rho[i] * e *(1-albedo[j]) * dE[j] ) / ( R[e] * W[e])
#    #%% Wedlund 2013 Eqn 1-3 (discretized form)
#    q = A.dot(Phi)
