#!/usr/bin/env python3
from __future__ import division,absolute_import
from datetime import datetime
from pytz import UTC
from numpy.testing import assert_allclose
#
from reesaurora.rees_model import reesiono,loadaltenergrid
from gridaurora.readApF107 import readmonthlyApF107


def test_reesiono():
    minalt=30; nalt=286
    isotropic=False
    glat=65; glon=-148
    mass=48.
    t = datetime(2013,3,31,12,tzinfo=UTC)
#%%
    apfn = 'data/RecentIndices.txt'
    f107Ap=readmonthlyApF107(int(str(t.year) + '{:02d}'.format(t.month)),
                             apfn)
    f107a = f107Ap['f107s']
    f107  = f107Ap['f107o']
    ap    = (f107Ap['Apo'],)*7
#%%
    z,E = loadaltenergrid(minalt,nalt)
#%%
    Q = reesiono(t,z,E,glat,glon,f107a,f107,ap,mass,isotropic)
#%%
    assert_allclose(Q.minor_axis.values,z)
    assert_allclose(Q.major_axis.values,E)
    assert Q.items.to_pydatetime()[0] == t

    Qv = Q.values
    assert_allclose([Qv[0,23,58],Qv[0,53,88]],[9.0995363363765226e-06,2.1067761704173373e-07])

if __name__ == '__main__':
    test_reesiono()