from matplotlib.pyplot import figure,subplots
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

def fig11(E,chi,Lambda_m,Lambda_i): #from Sergienko & Ivanov 1993
    # Lambda: Nenergy x Nalt

    def f11(Lambda,ax,ttxt):
        for e,x,l in zip(E,chi,Lambda):
            ax.plot(x,l,label=str(e))
        ax.set_title('$\lambda$ '+ttxt)
        ax.set_xlabel('Distance $\chi$')
        ax.legend()

    fg,axs = subplots(1,2,sharey=True)
    axs[0].set_ylabel('Function of Dissipation')

    f11(Lambda_m,axs[0],'field-aligned flux')

    f11(Lambda_i,axs[1],'isotropic flux')

