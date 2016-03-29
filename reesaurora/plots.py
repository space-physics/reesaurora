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

    def _f11(Lambda,ax,ttxt):
        for e,x,l in zip(E,chi,Lambda):
            ax.plot(x,l,label=str(e))
        ax.set_title('$\lambda$ '+ttxt)
        ax.set_xlabel('Distance $\chi$')
        ax.legend()

    fg,axs = subplots(1,2,sharey=True,num=11)
    axs[0].set_ylabel('Function of Dissipation')

    _f11(Lambda_m,axs[0],'field-aligned flux')

    _f11(Lambda_i,axs[1],'isotropic flux')

def fig12(E,Cmono,Ciso):
    """
    Energy deposition fitted curves
    """
    def f12(ax,C,E):
        for c,l in zip(C,('C1','C2','C3','C4')):
            ax.plot(E,c,label=l)
        ax.legend()
        ax.grid(which='both')
        ax.set_ylabel('fitting parameter')
        ax.set_xscale('log')
        ax.grid(True,which='both')
        ax.legend()
        ax.autoscale(True,tight=True, axis='x')

    fg,axs = subplots(2,1,sharex=True,num=12)

    ax=axs[0]
    f12(ax,Cmono,E)
    ax.set_title('Monodirectional')
#%%
    ax = axs[1]
    f12(ax,Ciso,E)
    ax.set_xlabel('Energy [eV]')
    ax.set_title('Isotropic Downward Hemisphere')

def fig13(E,afmono,afiso,rngmono,rngiso):

    fg,axs = subplots(2,1,sharex=True,num=13)
#%% range
    ax = axs[0]
    ax.plot(E,rngmono,label='monodirectional')
    ax.plot(E,rngiso,label='isotropic')
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Range [g/cm^2]')
    ax.legend(loc='best')
    ax.autoscale(True,tight=True,axis='x')
    ax.set_yscale('log')
    ax.set_ylim(1e-7,1e-4)
    ax.grid(which='both')
#%% albedo
    ax = axs[1]
    ax.plot(E,afmono,label='monodirectional')
    ax.plot(E,afiso,label='isotropic')
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Albedo-flux')
    ax.set_xscale('log')
    ax.legend(loc='best')
    ax.autoscale(True,tight=True,axis='x')