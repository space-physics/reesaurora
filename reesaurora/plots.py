from matplotlib.pyplot import figure, subplots
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator


def plotA(Q, ttxt, vlim):
    E = Q.energy.values
    z = Q.altkm.values
    # species = Q.species.values
    Q = Q[0]  # FIXME first time only
# %%

    def _doax(fg, ax, ttxt):
        ax.yaxis.set_major_locator(MultipleLocator(100))
        ax.yaxis.set_minor_locator(MultipleLocator(20))
        ax.set_xscale('log')
        ax.set_title(ttxt)
        ax.set_xlabel('beam energy [eV]')
        ax.autoscale(True, tight=True)  # fill axes

    # fg,axs = subplots(1,Q.species.size,sharey=True)
    # fg.suptitle(ttxt)
    fg = figure()

    # for ax,Qs,l in zip(axs,Q.values,species):
    Qs = Q.values
    ax = fg.gca()
    hi = ax.pcolormesh(E, z, Qs,
                       vmin=vlim[0], vmax=vlim[1],
                       norm=LogNorm())
    _doax(fg, ax, '')  # str(l))

    c = fg.colorbar(hi, ax=ax)
    c.set_label('Production Rate [cm$^{-3}$ s${^-1}$]')

    ax.set_ylabel('altitude [km]')


def fig7(W, z, Eplot):
    """
    Eqn. 6 uses each species individually, while Eqn A1 (actually used in the program)
    considers overall mass density of all three species at once.
    """
    # W: Nenergy x Nalt
    ax = figure().gca()
    for i, e in enumerate(Eplot):
        ax.plot(W[i, :], z, label='{:.0f}'.format(e))
    ax.set_title('Eqn A1')
    ax.set_xscale('log')
    ax.set_ylim(80, 275)
    ax.set_xlim(1e7, 1e9)
    ax.set_xlabel('Energy Deposition [eV g$^{-1}$ cm$^2$]')
    ax.set_ylabel('altitude [km]')
    ax.legend()


def fig8(Q):
    Q = Q.squeeze()
    E = Q.energy.values
    z = Q.altkm.values
    # species=Q.species.values

    # fg,axs = subplots(1,3,num=8,sharey=True)
    fg = figure()
    # for ax,s in zip(axs,species):
    ax = fg.gca()
    for e in E:
        ax.plot(Q.loc[:, e], z, label=str(e))
    ax.legend()
    ax.grid(True, which='both')
    ax.set_ylim(80, 275)
    ax.set_xscale('log')
    ax.set_xlim(1e-2, 1e5)
    ax.set_title('Ionization Rate')
    ax.set_xlabel('Ionization Rate [cm$^{-3}$ s$^{-1}$]')
    ax.set_ylabel('Altitude [km]')


def fig11(E, chi, Lambda_m, Lambda_i, fignum=11):  # from Sergienko & Ivanov 1993
    # Lambda: Nenergy x Nalt

    def _f11(Lambda, ax, ttxt):
        for e, x, l in zip(E, chi, Lambda):
            ax.plot(x, l, label=str(e))
        ax.set_title('$\lambda$ '+ttxt)
        ax.set_xlabel('Distance $\chi$')
        ax.set_xlim(0, 3)
        ax.legend()

    fg, axs = subplots(1, 2, sharey=True, num=fignum)
    axs[0].set_ylabel('Function of Dissipation')

    if Lambda_m is not None:
        _f11(Lambda_m, axs[0], 'field-aligned flux')

    if Lambda_i is not None:
        _f11(Lambda_i, axs[1], 'isotropic flux')


def fig12(E, Cmono, Ciso):
    """
    Energy deposition fitted curves
    """
    def f12(ax, C, E):
        for c, l in zip(C, ('C1', 'C2', 'C3', 'C4')):
            ax.plot(E, c, label=l)
        ax.legend()
        ax.grid(which='both')
        ax.set_ylabel('fitting parameter')
        ax.set_xscale('log')
        ax.grid(True, which='both')
        ax.legend()
        ax.autoscale(True, tight=True, axis='x')

    fg, axs = subplots(2, 1, sharex=True, num=12)

    ax = axs[0]
    f12(ax, Cmono, E)
    ax.set_title('Monodirectional')
# %%
    ax = axs[1]
    f12(ax, Ciso, E)
    ax.set_xlabel('Energy [eV]')
    ax.set_title('Isotropic Downward Hemisphere')


def fig13(E, afmono, afiso, rngmono, rngiso):

    fg, axs = subplots(2, 1, sharex=True, num=13)
# %% range
    ax = axs[0]
    ax.plot(E, rngmono, label='monodirectional')
    ax.plot(E, rngiso, label='isotropic')
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Range [g/cm^2]')
    ax.legend(loc='best')
    ax.autoscale(True, tight=True, axis='x')
    ax.set_yscale('log')
    ax.set_ylim(1e-7, 1e-4)
    ax.grid(which='both')
# %% albedo
    ax = axs[1]
    ax.plot(E, afmono, label='monodirectional')
    ax.plot(E, afiso, label='isotropic')
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Albedo-flux')
    ax.set_xscale('log')
    ax.legend(loc='best')
    ax.autoscale(True, tight=True, axis='x')
