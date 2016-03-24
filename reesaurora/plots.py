from matplotlib.pyplot import figure
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator

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

def fig11(E,Lambda):
    # Lambda: Nenergy x Nalt

    ax = figure().gca()
    for e,l in zip(E,Lambda):
        ax.plot(l,label=str(e))
    ax.set_title('$\lambda$')
    ax.legend()