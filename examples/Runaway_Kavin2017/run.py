#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
import sys
sys.path.append('../Deuterium_Burn-through')
sys.path.append('../../py')

#from run import makeplots
import PlasmaParameters as Formulas
import DREAM.Settings.Equations.ColdElectronTemperature as Tcold
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver

import DREAM.Settings.Atomics as Atomics
from STREAM import STREAMOutput, STREAMSettings, runiface
import STREAM.Settings.Equations.ElectricField as ElectricField
import STREAM.Settings.Equations.IonSpecies as Ions



def generate(prefill=3e-4, gamma=3e-2, Vloop=10.6, Vloop_t=0, Ures0=14, tmax=0.003, nt=1000):
    """
    Generate a STREAMSettings object for a simulation with the specified
    parameters.
    
    :param prefill: Prefill gas pressure [Torr]
    :param gamma:   Degree of ionization
    :param Vloop:   External loop voltage [V]
    :param Vloop_t: Time vector corresponding to given loop voltage [V]
    :param j0:      Initial plasma current density [A/m^2]
    :param tmax:    Simulation time [s]
    """

    '''
    #V_v     = 1700     # [m^3] (Vessel Volume)
    #R       = 5.7      # [m]   (major radius)
    #a       = 1.6      # [m]   (minor radius)
    #B       = 2.65     # [T]   (Toroidal magnetic field, at R = 6.2 m)
    #Y_H     = 1        #       (Recycling coefficient for hydrogen)
    #U_ext   = 10.6     # [V]   (External voltage)
    #T_e0    = 5        # [eV]  (Initial electron temperature)
    #T_i0    = 1        # [eV]  (Initial ion temperature)
    #gamma_0 = 0.03     #       (Initial degree of ionization)
    #t_bd    = 0.7      # [s]   (Time of gas breakdown, simulation starts)
    #t_end   = 2        # [s]   (Simulation ends)
    '''

    n0 = 0.01e20/7#3.22e22 * prefill  # Initial total deuterium density
    # n0 = 2.78e22 * prefill  # Initial total deuterium density
    nD = n0 * np.array([[1 - gamma], [gamma]])

    Btor = 2.65     # Toroidal magnetic field [T]
    a = 1.6         # Plasma minor radius [m]
    R0 = 5.7        # Plasma major radius [m]
    l_MK2 = 1       # Distance between plasma centre and passive structure [m] (unused)
    V_vessel = 1700 # Vacuum vessel volume [m^3]

    Te0 = 5     # electron temperature [eV]
    Ti0 = 1     # ion temperature [eV]

    # Initial electric field
    #E0 = j0 / Formulas.evaluateSpitzerConductivity(n=nD[1], T=Te0, Z=1)
    E0 = Ures0/(2 * np.pi * R0)

    ss = STREAMSettings()

    ss.atomic.adas_interpolation = Atomics.ADAS_INTERP_BILINEAR

    # Electric field
    ss.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
    ss.eqsys.E_field.setInitialProfile(E0)
    Lp = float(scipy.constants.mu_0 * R0 * (np.log(8*R0/a) + 0.5 - 2))
    ss.eqsys.E_field.setInductances(Lp=Lp, Lwall=9.1e-6, M=2.49e-6, Rwall=1e6)
    ss.eqsys.E_field.setCircuitVloop(Vloop, Vloop_t)

    # Electron temperature
    ss.eqsys.T_cold.setType(Tcold.TYPE_SELFCONSISTENT)
    ss.eqsys.T_cold.setInitialProfile(Te0)
    ss.eqsys.T_cold.setRecombinationRadiation(recombination=Tcold.RECOMBINATION_RADIATION_INCLUDED)

    # Ions
    ss.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=nD, r=np.array([0]), T=Ti0)
    ss.eqsys.n_i.addIon(name='Fe', Z=26, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=0, r=np.array([0]), T=Ti0)

    # Enable runaway
    ss.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_NEGLECT)
    ss.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)

    # Recycling coefficients (unused)
    ss.eqsys.n_i.setJET_CWrecycling()

    # Radial grid
    ss.radialgrid.setB0(Btor)
    ss.radialgrid.setMinorRadius(a)
    ss.radialgrid.setMajorRadius(R0)
    ss.radialgrid.setWallRadius(a)
    ss.radialgrid.setVesselVolume(V_vessel)

    ss.radialgrid.setRecyclingCoefficient1(1)
    ss.radialgrid.setRecyclingCoefficient2(0)
    ss.radialgrid.setRecyclingCoefficient3(1)
    
    # Disable kinetic grids
    ss.hottailgrid.setEnabled(False)
    ss.runawaygrid.setEnabled(False)

    # Numerical settings
    ss.solver.setType(Solver.NONLINEAR)
    ss.solver.preconditioner.setEnabled(False)
    ss.timestep.setTmax(tmax)
    ss.timestep.setNt(nt)
    ss.timestep.setNumberOfSaveSteps(10000)

    ss.other.include('fluid', 'stream', 'scalar')

    return ss


def drawplot1(axs, so, toffset=0.7):
    """
    Draw a plot with a output from the given STREAMOutput object.
    """
    t = so.grid.t[:] + toffset

    Te = so.eqsys.T_cold[:, 0]
    Ti = so.eqsys.W_i.getTemperature()['D'][:, 0]

    ne = so.eqsys.n_cold[:, 0]
    nD0 = so.eqsys.n_i['D'][0][:]
    nD1 = so.eqsys.n_i['D'][1][:]

    Vp = so.other.stream.V_p[:, 0] / 1e6
    Poh = -so.other.fluid.Tcold_ohmic[:, 0] * Vp
    Prad = so.other.fluid.Tcold_radiation[:, 0] * Vp

    Uext = 10.6 * t / t
    Ures = 2 * np.pi * 5.7 * so.eqsys.E_field[:,0]

    Ip = so.eqsys.I_p[:,0]
    Ire = 1.602e-19 * 2.99792458e8 * 1.6**2 * np.pi * so.eqsys.n_re[:]

    plotInternal(axs[0, 0], t, Te / 1e3, ylabel=r'$T$ (keV)', color='k', showlabel=True, label='$T_{\rm e}$')
    plotInternal(axs[0, 0], t, Ti / 1e3, ylabel=r'$T$ (keV)', color='m', showlabel=True, label='$T_{\rm i}$')

    plotInternal(axs[0, 1], t, ne / 1e20, ylabel=r'$n$ (1e20 m$^{-3}$)', color='k', showlabel=True, label='$n_{\rm e}$')
    plotInternal(axs[0, 1], t, 7 * nD0 / 1e20, ylabel=r'$n$ (1e20 m$^{-3}$)', color='m', showlabel=True, label='$n_{\rm D0}$')
    plotInternal(axs[0, 1], t, nD1 / 1e20, ylabel=r'$n$ (1e20 m$^{-3}$)', color='c', showlabel=True, label='$n_{\rm D1}$')

    plotInternal(axs[1, 0], t[1:], Poh, ylabel=r'$P$ (MW)', color='m', showlabel=True, label='$P_{\rm oh}$')
    plotInternal(axs[1, 0], t[1:], Prad, ylabel=r'$P$ (MW)', color='r', showlabel=True, label='$P_{\rm rad}$')

    plotInternal(axs[1, 1], t, Uext, ylabel=r'$U$ (V)', color='k', showlabel=True, label='$U_{\rm ext}$')
    plotInternal(axs[1, 1], t, Ures, ylabel=r'$U$ (V)', color='m', showlabel=True, label='$U_{\rm res}$')

    plotInternal(axs[2, 0], t, Ip / 1e6, ylabel=r'$I$ (MA)', color='k', showlabel=True, label='$I_{\rm p}$')
    plotInternal(axs[2, 0], t, Ire / 1e6, ylabel=r'$I$ (MA)', color='m', showlabel=True, label='$I_{\rm re}$')

    for i in range(axs.shape[0]):
        for j in range(axs.shape[1]):
            axs[i,j].set_xlim([t[0], t[-1]])
            axs[i,j].grid(True)
'''
    axs[0,0].set_ylim([0, 0.1])
    axs[0,1].set_ylim([0, 0.012])
    axs[1,0].set_ylim([0, 2.5])
    axs[1,1].set_ylim([0, 16])
    axs[2,0].set_ylim([0, 0.7])
'''
    #axs[0,0].set_yticks([0, 50, 100, 150, 200])
    #axs[0,1].set_yticks([0, 5, 10, 15])
    #axs[1,0].set_yticks([0, 20, 40, 60, 80])
    #axs[2,0].set_yticks([0, 20, 40, 60, 80])
    #axs[2,1].set_yticks([0, 0.1, 0.2, 0.3, 0.4])


def drawplot2(axs, so, toffset=0):
    t = so.grid.t[:] + toffset
    showlabel = (toffset==0)

    Vp = so.other.stream.V_p[:,0] / 1e6
    Poh = -so.other.fluid.Tcold_ohmic[:,0] * Vp
    Pequi = so.other.fluid.Tcold_ion_coll[:,0] * Vp
    Pconve = so.other.scalar.energyloss_T_cold[:,0] / 1e6
    Prad = so.other.fluid.Tcold_radiation[:,0] * Vp
    Pre = so.other.fluid.Tcold_nre_coll[:,0] * Vp

    PequiI = so.other.stream.Wi_e_coll[:,0] * Vp
    Pcx = -so.other.stream.Wi_chargeexchange[:,0] * Vp
    Pconvi = -so.other.stream.Wi_iontransport[:,0] * Vp

    ylbl = 'MW'
    plotInternal(axs[0], t[1:], Poh, ylabel=ylbl, color='r', linestyle='--', label='Ohmic', showlabel=showlabel)
    plotInternal(axs[0], t[1:], Pequi, ylabel=ylbl, color='b', linestyle='--', label='Equilibration', showlabel=showlabel)
    plotInternal(axs[0], t[1:], Pconve, ylabel=ylbl, color='k', linestyle='--', label='Transport', showlabel=showlabel)
    plotInternal(axs[0], t[1:], Prad, ylabel=ylbl, color='m', linestyle='--', label='Radiation + ionization', showlabel=showlabel)
    plotInternal(axs[0], t[1:], Pre, ylabel=ylbl, color='y', linestyle='--', label='Runaway collisions', showlabel=showlabel)

    plotInternal(axs[1], t[1:], PequiI, color='b', linestyle='--', ylabel=ylbl, label='Equilibration', showlabel=showlabel)
    plotInternal(axs[1], t[1:], Pcx, color='r', linestyle='--', ylabel=ylbl, label='Charge exchange', showlabel=showlabel)
    plotInternal(axs[1], t[1:], Pconvi, ylabel=ylbl, color='k', linestyle='--', label='Transport', showlabel=showlabel)

    axs[0].set_title('Electron energy balance')
    axs[1].set_title('Ion energy balance')
    for ax in axs:
        ax.set_xlim([0, 0.3])
        ax.set_ylim([0, 0.5])
        ax.legend(frameon=False)
        ax.grid(True)
        ax.set_xticks([0, 0.05, 0.1, 0.15, 0.2, 0.25])
        ax.set_yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5])


def drawplot3(axs, so, toffset=0):
    t = so.grid.t[:] + toffset
    showlabel = (toffset==0)

    axs[0,0].set_title('Plasma current')
    plotInternal(axs[0,0], t, so.eqsys.I_p[:,0]/1e3, ylabel='kA', color='k', linestyle='-', label='Total', showlabel=showlabel)
    plotInternal(axs[0,0], t, so.eqsys.j_re.current()/1e3, ylabel='kA', color='r', linestyle='--', label='Runaway', showlabel=showlabel)
    plotInternal(axs[0,0], t, so.eqsys.j_ohm.current()/1e3, ylabel='kA', color='b', linestyle=':', label='Ohmic', showlabel=showlabel)

    ylbl = r's$^{-1}$'
    rrNet = so.other.fluid.runawayRate[:,0]
    rrDreicer = so.other.fluid.gammaDreicer[:,0]
    rrAva = so.other.fluid.GammaAva[:,0]*so.eqsys.n_re[1:,0]
    rrTransp = rrDreicer + rrAva - rrNet

    axs[0,1].set_title('Runaway rate')
    plotInternal(axs[0,1], t[1:], rrNet, ylabel=ylbl, color='k', linestyle='-', label='Net', showlabel=showlabel)
    plotInternal(axs[0,1], t[1:], rrDreicer, ylabel=ylbl, color='b', linestyle=':', label='Dreicer', showlabel=showlabel)
    plotInternal(axs[0,1], t[1:], rrAva, ylabel=ylbl, color='r', linestyle='--', label='Avalanche', showlabel=showlabel)
    plotInternal(axs[0,1], t[1:], rrTransp, ylabel=ylbl, color='g', linestyle='-.', label='Avalanche', showlabel=showlabel)

    axs[1,0].set_title('Electric field')
    ylbl = r'$E_{\rm D}$ (\%)'
    plotInternal(axs[1,0], t[1:], so.eqsys.E_field.norm('ED')[1:,0] * 100, ylabel=ylbl, color='k', linestyle='-', label=r'$E/E_{\rm D}$', showlabel=showlabel)
    plotInternal(axs[1,0], t[1:], so.other.fluid.Eceff[:,0] / so.other.fluid.EDreic[:,0] * 100, ylabel=ylbl, color='r', linestyle='--', label=r'$E_{\rm c,eff}$', showlabel=showlabel)
    axs[1,0].set_ylim([0, 10])

    axs[1,1].set_title('Confinement time')
    ylbl = r's'
    plotInternal(axs[1,1], t[1:], so.other.stream.tau_RE[:,0], ylabel=ylbl, color='k', linestyle='-', label=r'$\tau_{\rm RE}$', showlabel=showlabel)
    plotInternal(axs[1,1], t[1:], so.other.stream.tau_RE1[:,0], ylabel=ylbl, color='r', linestyle='--', label=r'$\tau_{\rm RE,1}$', showlabel=showlabel)
    plotInternal(axs[1,1], t[1:], so.other.stream.tau_RE2[:,0], ylabel=ylbl, color='b', linestyle=':', label=r'$\tau_{\rm RE,2}$', showlabel=showlabel)
    axs[1,1].set_ylim([0, 0.1])

    for ax in axs.flatten():
        ax.set_xlim([0, t[-1]])
        ax.legend(frameon=False)
        ax.grid(True)


def plotInternal(ax, x, y, ylabel, xlbl=True, ylim=None, log=False, showlabel=False, label=None, *args, **kwargs):
    if label is not None and showlabel == False:
        label = None

    if log:
        ax.semilogy(x, y, label=label, *args, **kwargs)
    else:
        ax.plot(x, y, label=label, *args, **kwargs)

    ax.set_xlim([0, x[-1]])
    if xlbl:
        ax.set_xlabel(r'Time (s)')
    ax.set_ylabel(ylabel)

    if ylim is not None:
        ax.set_ylim(ylim)


def makeplots(so1, so2):
    fig1, axs1 = plt.subplots(3, 2, figsize=(7, 10))

    drawplot1(axs1, so1)
    drawplot1(axs1, so2, toffset=so1.grid.t[-1])

    #fig2, axs2 = plt.subplots(1, 2, figsize=(10, 4))

    #drawplot2(axs2, so1)
    #drawplot2(axs2, so2, toffset=so1.grid.t[-1])

    #fig3, axs3 = plt.subplots(2, 2, figsize=(10, 8))

    #drawplot3(axs3, so1)
    #drawplot3(axs3, so2, toffset=so1.grid.t[-1])

    fig1.tight_layout()
    #fig2.tight_layout()
    #fig3.tight_layout()
    plt.show()


def main(argv):
    FONTSIZE = 16
    plt.rcParams.update({'font.size': FONTSIZE})

    parser = argparse.ArgumentParser()

    parser.add_argument('-e', '--extension', help="Append the specified string to the end of file names", dest="extension", action='store', default='')
    parser.add_argument('-n', '--no-plot', help="Only run simulations, don't generate plots", dest="plot", action='store_false', default=True)
    parser.add_argument('-s', '--skip', help="Skip the simulation and load output from already existing files", action='store_true', dest="skip", default=False)

    settings = parser.parse_args()

    ext = '' if not settings.extension else '_' + settings.extension

    if not settings.skip:
        prefill = 3e-4
        ss1 = generate(prefill=prefill, nt=70000)
        ss1.save(f'settings1{ext}.h5')
        so1 = runiface(ss1, f'output1{ext}.h5', quiet=False)

        ss2 = STREAMSettings(ss1)
        ss2.fromOutput(f'output1{ext}.h5')
        ss2.timestep.setTmax(1.3 - ss1.timestep.tmax)
        ss2.timestep.setNumberOfSaveSteps(0)
        ss2.timestep.setNt(100000)
        ss2.save(f'settings2{ext}.h5')
        so2 = runiface(ss2, f'output2{ext}.h5', quiet=False)
    else:
        so1 = STREAMOutput(f'output1{ext}.h5')
        so2 = STREAMOutput(f'output2{ext}.h5')

    if settings.plot:
        makeplots(so1, so2)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


