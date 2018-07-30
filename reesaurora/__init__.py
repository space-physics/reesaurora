"""
 ionization_profiles_from_flux - simple model for volume emission as function of altitude.
   After Sergienko and Ivanov 1993
   a massively speeded up implementation after the AIDA_TOOLS package by Gustavsson, Brandstrom, et al
"""
from pathlib import Path
import logging
from dateutil.parser import parse
from datetime import datetime
import xarray
import numpy as np
from scipy.interpolate import interp1d  # not numpy.interp since only for 1-D
from typing import Union, Tuple
from msise00 import rungtd1d
from gridaurora.ztanh import setupz
from gridaurora.zglow import glowalt

species = ['N2', 'O', 'O2']
usesemeter = True


def reesiono(T: Union[str, datetime], altkm: np.ndarray, E: np.ndarray,
             glat: float, glon: float, isotropic: bool, verbose: bool,
             datfn: Path) -> xarray.DataArray:
    # other assertions covered inside modules
    assert isinstance(isotropic, bool)

    if abs(glat) < 45.:
        logging.error('This model was intended for auroral precipitation.')

    if isinstance(T, str):
        T = parse(T)
    time = np.atleast_1d(T)
    assert isinstance(time[0], datetime)
# %% MSIS
    if isotropic:
        logging.debug('isotropic pitch angle flux')
    else:
        logging.debug('field-aligned pitch angle flux')

    Qt = xarray.DataArray(data=np.empty((time.size, altkm.size, E.size)),
                          coords=[time, altkm, E],
                          dims=['time', 'alt_km', 'energy'])
# %% loop
    for t in time:
        iono = rungtd1d(t, altkm, glat, glon)

        Q = ionization_profile_from_flux(E, iono, isotropic, datfn, verbose)
        Qt.loc[t, ...] = Q

    return Qt


def ionization_profile_from_flux(E: np.ndarray, iono: xarray.Dataset,
                                 isotropic: bool, datfn: Path, verbose: bool) -> np.ndarray:
    """
    simple model for volume emission as function of altitude.
    After Sergienko and Ivanov 1993 and Gustavsson AIDA_TOOLs
    """
    if ((E < 50.) | (E > 1e4)).any():
        logging.warning('Sergienko & Ivanov 1993 covered E \in [100,10000] eV')

    if (iono.alt_km > 700.).any():
        logging.error('Sergienko & Ivanov 1993 assumed electron source was at altitude 700km.')

# %% Table 1 Sergienko & Ivanov 1993, rightmost column
    # mean energy per ion-electron pair
    E_cost_ion = {'N2': 36.8, 'O2': 28.2, 'O': 26.8}
# %% Eqn 7, Figure 6
    k = {'N2': 1., 'O2': 0.7, 'O': 0.4}

    dE = np.diff(E)
    dE = np.append(dE, dE[-1])

    Peps = partition(iono, k, E_cost_ion)  # Nalt x Nspecies

# %% First calculate the energy deposition as a function of altitude
    qE = np.empty((iono.alt_km.size, E.size))  # Nalt x Nenergy
    for i, (e, d) in enumerate(zip(E, dE)):
        Ebins = np.linspace(e, e+d, 20)
        # for isotropic or field aligned electron beams
        Am = energy_deg(Ebins, isotropic, iono)

        q = Am.sum(axis=0)  # sum over the interim energy sub-bins
        q *= Peps.sum(axis=1)  # effect of ion chemistry at each altitude
        qE[:, i] = q

    return qE


def energy_deg(E: np.ndarray, isotropic: bool, iono: xarray.Dataset) -> np.ndarray:
    """
    energy degradation of precipitating electrons
    """
    atmp = iono['Total'].squeeze() / 1e3

    N_alt0 = atmp.alt_km.size
    zetm = np.zeros(N_alt0)
    dH = np.gradient(atmp.alt_km)
    for i in range(N_alt0-1, 0, -1):  # careful with these indices!
        dzetm = (atmp[i]+atmp[i-1])*dH[i-1]*1e5/2
        zetm[i-1] = zetm[i] + dzetm

    alb = albedo(E, isotropic)
    assert E.shape == alb.shape

    dE = np.gradient(E)
    r = PitchAngle_range(E, isotropic)
    assert E.shape == r.shape

    hi = zetm / r[:, None]
    assert hi.shape == (E.size, zetm.size)

    Lambda = lambda_comp(hi, E, isotropic)

    Am = atmp.values * Lambda * E[:, None] * (1-alb[:, None]) / r[:, None]
    assert Am.shape == (E.size, zetm.size)

    Am[0, :] *= dE[0]/2.
    Am[-1, :] *= dE[-1]/2.
    Am[1:-2, :] *= (dE[1:-2] + dE[0:-3])[:, None]/2.

    return Am


def PitchAngle_range(E: np.ndarray, isotropic: bool) -> np.ndarray:
    pr = 1.64e-6 if isotropic else 2.16e-6
    return pr * (E/1e3)**1.67 * (1 + 9.48e-2 * E**-1.57)


def albedo(E: np.ndarray, isotropic: Union[int, bool]) -> np.ndarray:
    """ ionospheric albedo model"""
    isotropic = int(isotropic)
    assert isotropic in (0, 1)

    logE_p = np.append(1.69, np.arange(1.8, 3.7+0.1, 0.1))
    Param = np.array(
        [[0.352, 0.344, 0.334, 0.320, 0.300, 0.280, 0.260, 0.238, 0.218, 0.198,
          0.180, 0.160, 0.143, 0.127, 0.119, 0.113, 0.108, 0.104, 0.102, 0.101, 0.100],
         [0.500, 0.492, 0.484, 0.473, 0.463, 0.453, 0.443, 0.433, 0.423, 0.413,
          0.403, 0.395, 0.388, 0.379, 0.378, 0.377, 0.377, 0.377, 0.377, 0.377, 0.377]])
    logE = np.log10(E)

    falb = interp1d(logE_p, Param[isotropic, :], kind='linear', bounds_error=False, fill_value=np.nan)
    alb = falb(logE)
    alb[logE > logE_p[-1]] = Param[isotropic, -1]

    return alb


def lambda_comp(hi: np.ndarray, E: np.ndarray, isotropic: bool) -> np.ndarray:
    """
    interpolated over energies from 48.9 eV to 5012 eV
    for isotropic and field-aligned precipitation
    """

# %% field-aligned
    logE_m = np.append(1.69, np.arange(1.8, 3.7+0.1, 0.1))
    Param_m = np.array(
        [[1.43, 1.51, 1.58, 1.62, 1.51, 1.54, 1.18, 1.02, 0.85, 0.69, 0.52, 0.35,
          0.21, 0.104, 0.065, 0.05, 0.04, 0.03, 0.03, 0.025, 0.021],
         [0.83, 0.77, 0.72, 0.67, 0.63, 0.59, 0.56, 0.525, 0.495, 0.465, 0.44,
             0.42, 0.40, 0.386, 0.37, 0.36, 0.35, 0.34, 0.335, 0.325, 0.32],
         [-0.025, -0.030, -0.040, -0.067, -0.105, -0.155, -0.210, -0.275, -0.36, -0.445, -
             0.51, -0.61, -0.69, -0.77, -0.83, -0.865, -0.90, -0.92, -0.935, -0.958, -0.96],
         [-1.67, -1.65, -1.62, -1.56, -1.46, -1.35, -1.20, -0.98, -0.70, -0.37,
          -0.063, 0.39, 0.62, 0.92, 1.11, 1.25, 1.36, 1.44, 1.50, 1.55, 1.56]]
    )
# %% isotropic
    """
        interpolated over energies from 48.9 eV to 1000 eV
    """
    logE_i = np.append(1.69, np.arange(1.8, 3.0+0.1, 0.1))
    Param_i = np.array(
        [[0.041, 0.051, 0.0615, 0.071, 0.081, 0.09, 0.099, 0.1075, 0.116, 0.113, 0.13, 0.136, 0.139, 0.142],
         [1.07, 1.01, 0.965, 0.9, 0.845, 0.805, 0.77, 0.735, 0.71, 0.69, 0.67, 0.665, 0.66, 0.657],
         [-0.064, -0.1, -0.132, -0.171, -0.2, -0.221, -0.238, -0.252, -0.261, -0.267, -0.271, -0.274, -0.276, -0.277],
         [-1.054, -0.95, -0.845, -0.72, -0.63, -0.54, -0.475, -0.425, -0.38, -0.345, -0.319, -0.295, -0.28, -0.268]]
    )

    logE = np.log10(E)

    if isotropic:
        P = Param_i
        LE = logE_i
        Emax = 1000.
    else:
        P = Param_m
        LE = logE_m
        Emax = 5000.
# %% interpolate
    fC = interp1d(LE, P, kind='linear', axis=1, bounds_error=False, fill_value=np.nan)
    C = fC(logE)
# %% low energy
    lam = ((C[0, :][:, None]*hi + C[1, :][:, None]) *
           np.exp(C[2, :][:, None]*hi**2 + C[3, :][:, None]*hi))

    assert lam.shape == hi.shape
# %% high energy
    badind = E > Emax
    lam[badind] = (
        (P[0, -1]*hi[badind] + P[1, -1]) *
        np.exp(P[2, -1]*hi[badind]**2 + P[3, -1]*hi[badind])
    )

    return lam


def partition(iono: xarray.Dataset, k: np.ndarray, cost: np.ndarray) -> xarray.DataArray:
    """
    Implement Eqn 7 Sergienko 1993

    N: NUMBER density [cm^-3] vs. altitude for each species
    k: correction factors vs. Monte Carlo for Sergienko 1993
    cost: energization cost

    output:
    P_i(h) / epsilon
    """
    # m^-3 /1e6 = cm^-3
    N = iono[species] / 1e6  # [cm^-3]

    num = xarray.DataArray(data=np.empty((N.alt_km.size, len(species))),
                           coords=[N.alt_km, species],
                           dims=['alt_km', 'species'])
    for i in species:
        num.loc[:, i] = k[i] * N[i].squeeze()

    den = num.sum('species')

#    P = num.copy()
#    for i in species:
#        P.loc[:, i] = num.loc[:, i] / den
# %% energization cost
    Peps = num.copy()
    for i in species:
        Peps.loc[:, i] = num.loc[:, i] / den / cost[i]

    return Peps


def loadaltenergrid(minalt: float=90, Nalt: int=286, special_grid: str='') -> Tuple[np.ndarray, np.ndarray]:
    """
    makes a tanh-spaced grid (see setupz for info)

    minalt: [km] minimum altiude in grid (e.g. 90)
    Nalt: number of points in grid
    special_grid: use same grid as 'transcar' or 'glow'
    """
    assert isinstance(special_grid, str)
    # %% altitude
    if special_grid.lower() == 'transcar':
        z = setupz(286, 90, 1.5, 11.1475)
    elif special_grid.lower() == 'glow':
        z = glowalt()
    else:
        z = setupz(Nalt, minalt, 1.5, 11.1475)

    z = z[z <= 700]  # keeps original spacing, but with heights less than source at 700km
# %% energy of beams
    if special_grid.lower() == 'transcar':
        E = np.logspace(1.72, 4.25, num=33, base=10)
    else:
        E = np.logspace(1.72, 4.25, num=81, base=10)

    return z, E
