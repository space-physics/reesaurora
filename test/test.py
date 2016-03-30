#!/usr/bin/env python3
from pathlib import Path
from datetime import datetime
from numpy import datetime64
from numpy.testing import assert_allclose,run_module_suite
#
from reesaurora.rees_model import reesiono,loadaltenergrid

tdir  = Path(__file__).parents[1]

def test_reesiono():
    minalt=90; nalt=286
    isotropic=False
    glat=65; glon=-148
    t = '2013-03-31T12:00:00Z' #datetime(2013,3,31,12,tzinfo=UTC)
#%%
    z,E = loadaltenergrid(minalt,nalt) #altitude grid, Energy Grid
#%%
    Q = reesiono(t,z,E,glat,glon,isotropic,datfn=tdir/'data/SergienkoIvanov.h5')
#%%
    assert_allclose(Q.altkm.values,z)
    assert_allclose(Q.energy.values,E)
    assert Q.time[0] == datetime64(t)

    Qv = Q[0].sum(axis=0)
    assert_allclose([Qv[23,58],Qv[53,68]],
                    [0.0007949280872816374,4.2094326282647114e-06])

if __name__ == '__main__':
    run_module_suite()
