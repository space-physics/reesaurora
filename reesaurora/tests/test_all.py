from pathlib import Path
import xarray
from numpy import datetime64
from pytest import approx

import reesaurora as ra

R = Path(__file__).resolve().parents[1]


def test_reesiono():
    minalt = 90
    nalt = 286
    isotropic = False
    glat = 65.0
    glon = -148.0
    t = "2013-03-31T12:00:00"  # datetime(2013,3,31,12)
    # %%
    z, E = ra.loadaltenergrid(minalt, nalt)  # altitude grid, Energy Grid
    # %%
    Q = ra.reesiono(
        t,
        z,
        E,
        glat,
        glon,
        isotropic,
        datfn=R / "data/SergienkoIvanov.h5",
        verbose=False,
    )
    assert isinstance(Q, xarray.DataArray)
    # %%
    assert Q.alt_km.values == approx(z)
    assert Q.energy.values == approx(E)
    assert Q.time[0] == datetime64(t), "times didnt match up"

    Qv = Q.squeeze()
    # Qv = Q[0].sum('species')  # total production
    # print([Qv[23,58],Qv[53,68]])
    assert [Qv[23, 58], Qv[53, 68]] == approx([8.186955e-04, 5.609914e-06], rel=1e-5)
