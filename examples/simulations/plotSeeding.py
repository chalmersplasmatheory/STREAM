#!/usr/bin/env python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../py')

from STREAM import STREAMOutput


FONTSIZE = 20
mpl.rcParams.update({'text.usetex': True, 'font.family': 'sans', 'font.size': FONTSIZE})


def confgeneral(ax):
    ax.set_xlim([0, 8])
    ax.tick_params(axis='both', labelsize=FONTSIZE)


def confplotCurrent(ax):
    ax.set_ylim([0, 15])
    ax.set_xlabel(r'Time $t$ (s)', usetex=True)
    ax.set_ylabel(r'Plasma current $I_{\rm p}$ (MA)', usetex=True)
    ax.legend(frameon=False, fontsize=FONTSIZE)
    confgeneral(ax)


def confplotTemperature(ax):
    ax.set_ylim([0, 1.8e3])
    ax.set_xlabel(r'Time $t$ (s)', usetex=True)
    ax.set_ylabel(r'Electron temperature $T_{\rm e}$ (eV)', usetex=True)
    ax.set_yticks([0, 500, 1000, 1500])
    #ax.legend(frameon=False)
    confgeneral(ax)


def confplotEfield(ax):
    ax.set_ylim([0, 10])
    ax.set_xlabel(r'Time $t$ (s)', usetex=True)
    ax.set_ylabel(r'Electric field $E/E_{\rm D}$ (\%)', usetex=True)
    #ax.legend(frameon=False)
    confgeneral(ax)


def confplotRates(ax):
    ax.set_ylim([0, 10e15])
    ax.set_xlabel(r'Time $t$ (s)', usetex=True)
    ax.set_ylabel(r'$\partial n_{\rm re}/\partial t$ (m$^{-3}$s$^{-1}$)', usetex=True)
    ax.legend(frameon=False, fontsize=FONTSIZE)
    confgeneral(ax)


def confplotFuelling(ax):
    ax.set_ylim([0, 7e17])
    ax.set_xlabel(r'Time $t$ (s)', usetex=True)
    ax.set_ylabel(r'Fuelling rate $S_{D,{\rm fuel}}$ (m$^{-3}$s$^{-1}$)', usetex=True)
    confgeneral(ax)


def confplotDensity(ax):
    ax.set_ylim([0, 1.4*1e18])
    ax.set_xlabel(r'Time $t$ (s)', usetex=True)
    ax.set_ylabel(r'Electron density $n_{\rm e}$ (m$^{-3}$)', usetex=True)
    confgeneral(ax)


def plotCurrent(ax, so, color='k', label1=None, label2=None):
    h1, = ax.plot(so.grid.t, so.eqsys.I_p[:]/1e6, color=color, linestyle='-', label=label1)
    h2, = ax.plot(so.grid.t, so.eqsys.j_re.current()/1e6, color=color, linestyle='--', label=label2)

    return h1, h2


def plotTemperature(ax, so, color='k', label=None):
    ax.plot(so.grid.t, so.eqsys.T_cold[:], color=color, linestyle='-', label=label)
    #ax.plot(so.grid.t, so.eqsys.W_i.getTemperature('D')[:], color=color, linestyle='--')


def plotElectricField(ax, so, color='k', label=None):
    ax.plot(so.grid.t[1:], 100 * so.eqsys.E_field[1:,0] / so.other.fluid.EDreic[:,0], color=color, linestyle='-', label=label)


def plotDensity(ax, so, color='k', label=None):
    ax.plot(so.grid.t, so.eqsys.n_cold[:], color=color, linestyle='-', label=label)


def plotFuelling(ax, so, t1, t2, color='k'):
    t = so.grid.t[:]
    f = np.zeros(t.shape)
    f[np.where((t >= t1) & (t <= t2))] = 5e17

    h, = ax.plot(so.grid.t, f, color=color)

    return h


def plotRates(ax, so, color='k', label1=None, label2=None, label3=None):
    rr = so.other.fluid.runawayRate[:,0]
    gammaDreicer = so.other.fluid.gammaDreicer[:,0]
    GammaAva = so.other.fluid.GammaAva[:,0]
    nre = so.eqsys.n_re[1:,0]

    ax.plot(so.grid.t[1:], rr, color=color, linestyle='-', label=label1)
    ax.plot(so.grid.t[1:], gammaDreicer, color=color, linestyle='--', label=label2)
    ax.plot(so.grid.t[1:], GammaAva*nre, color=color, linestyle=':', label=label3)


def plotShading(ax, colors=['k', 'tab:red', 'tab:blue']):
    ax.axvspan(0.5, 2.5, color=colors[0], alpha=0.12, linewidth=0)
    ax.axvspan(1.0, 3.0, color=colors[1], alpha=0.12, linewidth=0)
    ax.axvspan(3.0, 5.0, color=colors[2], alpha=0.12, linewidth=0)


# Load simulation data
so05 = STREAMOutput('output/DwRE-fuel050-2.h5')
so10 = STREAMOutput('output/DwRE-fuel100-2.h5')
so30 = STREAMOutput('output/DwRE-fuel300-2.h5')

###########################
# FIRST FIGURE ("RUNAWAY")
###########################

fig1, axs1 = plt.subplots(1, 3, figsize=(14,4.5))

# Current
title = plotFuelling(axs1[0], so05, -2, -1, color='white')
h05, _ = plotCurrent(axs1[0], so05, color='k')
h10, _ = plotCurrent(axs1[0], so10, color='tab:red')
h30, _ = plotCurrent(axs1[0], so30, color='tab:blue')
plotCurrent(axs1[0], so05, color='k', label1=r'$I_{\rm p}$', label2=r'$I_{\rm re}$')
plotShading(axs1[0])
confplotCurrent(axs1[0])

# Electric field
plotElectricField(axs1[1], so05, color='k')
plotElectricField(axs1[1], so10, color='tab:red')
plotElectricField(axs1[1], so30, color='tab:blue')
plotShading(axs1[1])
confplotEfield(axs1[1])

# Runaway rates
plotRates(axs1[2], so05, color='k', label1='Total', label2='Dreicer', label3='Avalanche')
plotRates(axs1[2], so10, color='tab:red')
plotRates(axs1[2], so30, color='tab:blue')
plotShading(axs1[2])
confplotRates(axs1[2])

fig1.legend([title, h05, h10, h30], ['Fuelling duration:', '$0.5$-$2.5\,\mathrm{s}$', '$1$-$3\,\mathrm{s}$', '$3$-$5\,\mathrm{s}$'], loc='upper center', frameon=False, ncol=4)
fig1.tight_layout()
fig1.subplots_adjust(top=0.83)

############################
# SECOND FIGURE ("PLASMA")
############################
fig2, axs2 = plt.subplots(1, 3, figsize=(14,4.5))

# Fuelling function
title = plotFuelling(axs2[0], so05, -2, -1, color='white')
h05 = plotFuelling(axs2[0], so05, 0.5, 2.5, color='k')
h10 = plotFuelling(axs2[0], so10, 1.0, 3.0, color='tab:red')
h30 = plotFuelling(axs2[0], so30, 3.0, 5.0, color='tab:blue')
confplotFuelling(axs2[0])

plotDensity(axs2[1], so05, color='k')
plotDensity(axs2[1], so10, color='tab:red')
plotDensity(axs2[1], so30, color='tab:blue')
plotShading(axs2[1])
confplotDensity(axs2[1])

# Temperature
plotTemperature(axs2[2], so05, color='k')
plotTemperature(axs2[2], so10, color='tab:red')
plotTemperature(axs2[2], so30, color='tab:blue')
plotShading(axs2[2])
confplotTemperature(axs2[2])

fig2.legend([title, h05, h10, h30], ['Fuelling duration:', '$0.5$-$2.5\,\mathrm{s}$', '$1$-$3\,\mathrm{s}$', '$3$-$5\,\mathrm{s}$'], loc='upper center', frameon=False, ncol=4)
fig2.tight_layout()
fig2.subplots_adjust(top=0.83)

plt.show()

