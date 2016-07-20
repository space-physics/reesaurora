#!/usr/bin/env python
from reesaurora import Path
from numpy import datetime64
from numpy.testing import assert_allclose,run_module_suite
#
from reesaurora.rees_model import reesiono,loadaltenergrid

tdir  = Path(__file__).parents[1]

def test_reesiono():
    minalt=90; nalt=286
    isotropic=False
    glat=65.; glon=-148.
    t = '2013-03-31T12:00:00Z' #datetime(2013,3,31,12,tzinfo=UTC)
#%%
    z,E = loadaltenergrid(minalt,nalt) #altitude grid, Energy Grid
#%%
    Q = reesiono(t,z,E,glat,glon,isotropic,datfn=tdir/'data/SergienkoIvanov.h5',verbose=0)
#%%
    assert_allclose(Q.altkm.values,z)
    assert_allclose(Q.energy.values,E)
    assert Q.time[0] == datetime64(t),'times didnt match up'

    Qv = Q.squeeze()
    #Qv = Q[0].sum('species')  # total production
    #print([Qv[23,58],Qv[53,68]])
    assert_allclose([Qv[23,58],Qv[53,68]],
                    [8.280217e-04,   5.448110e-06],rtol=1e-5)

if __name__ == '__main__':
    run_module_suite()
