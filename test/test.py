#!/usr/bin/env python3
from datetime import datetime
from pytz import UTC
from numpy.testing import assert_allclose,run_module_suite
#
from reesaurora.rees_model import reesiono,loadaltenergrid

def test_reesiono():
    minalt=30; nalt=286
    isotropic=False
    glat=65; glon=-148
    t = '2013-03-31T12:00:00Z' #datetime(2013,3,31,12,tzinfo=UTC)
#%%
    z,E = loadaltenergrid(minalt,nalt)
#%%
    Q = reesiono(t,z,E,glat,glon,isotropic)
#%%
    assert_allclose(Q.minor_axis.values,z)
    assert_allclose(Q.major_axis.values,E)
    assert Q.items.to_pydatetime()[0] == datetime(2013,3,31,12,tzinfo=UTC)

    Qv = Q.values
    assert_allclose([Qv[0,23,58],Qv[0,53,88]],[9.0995363363765226e-06,2.1067761704173373e-07])

if __name__ == '__main__':
    run_module_suite()