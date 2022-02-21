#!/usr/bin/env python3
#
# Script to reproduce the ITER simulations of section 2 in (Kim et al., NF 2020).
# Note that in order for agreement between STREAM and DYON/BKD0/SCENPLINT to
# be achieved, it is necessary to use C-X rates for deuterium, even though the
# main ion species used here is hydrogen. This can only be achieved by modifying
# the source code as has been done on the 'Kim2020ITER' branch of STREAM.

import argparse
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
import sys
sys.path.append('../../Deuterium_Burn-through')
sys.path.append('../../../py')

#from run import makeplots
import PlasmaParameters as Formulas
import DREAM.Settings.Equations.ColdElectronTemperature as Tcold
import DREAM.Settings.Solver as Solver
import DREAM.Settings.Equations.OhmicCurrent as JOhm

import DREAM.Settings.Atomics as Atomics
from STREAM import STREAMOutput, STREAMSettings, runiface
import STREAM.Settings.Equations.ElectricField as ElectricField
import STREAM.Settings.Equations.IonSpecies as Ions


def generate(prefill=5e-5, gamma=2e-3, Vloop=12, Vloop_t=0, j0=298.4, tmax=0.003, nt=1000):
    """
    Generate a STREAMSettings object for a simulation with the specified
    parameters.
    
    :param prefill: Prefill gas pressure [Torr]
    :param gamma:   Degree of ionization
    :param Vloop:   External loop voltage [V]
    :param Vloop_t: Time vector corresponding to given loop voltage [V]
    :param j0:      Initial plasma current density [A/m^2]
    :param tmax:    Simulation time [s]
    :param nt:      Number of time steps
    """
    n0 = 3.22e22 * prefill  # Initial total deuterium density
    #n0 = 2.78e22 * prefill  # Initial total deuterium density
    nD = n0 * np.array([[1-gamma], [gamma]])

    a = 1.6         # Plasma minor radius [m]
    R0 = 5.65       # Plasma major radius [m]
    Btor = 2.65 * 6.2 / R0  # Toroidal magnetic field [T]
    #Btor = 2.65
    l_MK2 = 1       # Distance between plasma centre and passive structure [m] (unused)
    V_vessel = 1000 # Vacuum vessel volume

    Te0 = 1     # electron temperature [eV]
    Ti0 = 0.026 # ion temperature [eV]

    # Initial electric field
    E0 = j0 / Formulas.evaluateSpitzerConductivity(n=nD[1], T=Te0, Z=1)

    ss = STREAMSettings()

    ss.atomic.adas_interpolation = Atomics.ADAS_INTERP_BILINEAR

    # Electric field
    ss.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
    ss.eqsys.E_field.setInitialProfile(E0)
    #ss.eqsys.E_field.setInductances(Lp=6.09e-6, Lwall=9.1e-6, M=2.49e-6, Rwall=1e6)
    Lp = float(scipy.constants.mu_0 * R0 * (np.log(8*R0/a) + 0.25 - 2))
    ss.eqsys.E_field.setInductances(Lp=Lp, Lwall=9.1e-6, M=2.49e-6, Rwall=1e6)
    ss.eqsys.E_field.setCircuitVloop(Vloop, Vloop_t)
    #ss.eqsys.j_ohm.setConductivityMode(JOhm.CONDUCTIVITY_MODE_BRAAMS)
    ss.eqsys.j_ohm.setConductivityMode(JOhm.CONDUCTIVITY_MODE_SAUTER_COLLISIONAL)

    # Electron temperature
    ss.eqsys.T_cold.setType(Tcold.TYPE_SELFCONSISTENT)
    ss.eqsys.T_cold.setInitialProfile(Te0)
    ss.eqsys.T_cold.setRecombinationRadiation(recombination=Tcold.RECOMBINATION_RADIATION_INCLUDED)

    # Ions
    ss.eqsys.n_i.addIon(name='H', Z=1, iontype=Ions.IONS_DYNAMIC, n=nD, r=np.array([0]), T=Ti0, hydrogen=True)

    # Disable runaway
    ss.eqsys.n_re.setAvalanche(False)
    ss.eqsys.n_re.setDreicer(False)

    # Recycling coefficients (unused)
    ss.eqsys.n_i.setJET_CWrecycling()

    # Radial grid
    ss.radialgrid.setB0(Btor)
    ss.radialgrid.setMinorRadius(a)
    ss.radialgrid.setMajorRadius(R0)
    ss.radialgrid.setWallRadius(a)
    ss.radialgrid.setVesselVolume(V_vessel)
    ss.radialgrid.setBv(2e-3)

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


def drawplot1(axs, so, toffset=0, showlabel=False, save=True):
    """
    Draw a plot with a output from the given STREAMOutput object.
    """
    t = so.grid.t[:] + toffset
    Ip = so.eqsys.I_p[:, 0]
    ne = so.eqsys.n_cold[:, 0]
    Te = so.eqsys.T_cold[:, 0]
    #Ti = so.eqsys.W_i.getTemperature()['D'][:, 0]
    Ti = so.eqsys.W_i.getTemperature()['H'][:, 0]
    Lf = so.other.stream.Lf[:, 0]
    tau = so.other.stream.tau_D[:, 0]

    if save:
        t_csv = open('Data/time_STREAM.csv', 'ab')
        np.savetxt(t_csv, t)
        t_csv.close()
        Ip_csv = open('Data/PlasmaCurrent_STREAM.csv', 'ab')
        np.savetxt(Ip_csv, Ip)
        Ip_csv.close()
        ne_csv = open('Data/ElectronDensity_STREAM.csv', 'ab')
        np.savetxt(ne_csv, ne)
        ne_csv.close()
        Te_csv = open('Data/ElectronTemperature_STREAM.csv', 'ab')
        np.savetxt(Te_csv, Te)
        Te_csv.close()
        Ti_csv = open('Data/IonTemperature_STREAM.csv', 'ab')
        np.savetxt(Ti_csv, Ti)
        Ti_csv.close()
        Lf_csv = open('Data/ConnectionLength_STREAM.csv', 'ab')
        np.savetxt(Lf_csv, Lf)
        Lf_csv.close()
        tau_csv = open('Data/ConfinementTime_STREAM.csv', 'ab')
        np.savetxt(tau_csv, tau)
        tau_csv.close()

    plotInternal(axs[0,0], t, Ip/1e3, ylabel=r'$I_{\rm p}$ (kA)', color='g', showlabel=showlabel, label='STREAM')
    plotInternal(axs[0,1], t, ne/1e17, ylabel=r'$n_{\rm e}$ ($1\cdot 10^{17}$m$^{-3}$)', color='g', showlabel=False, label='STREAM')
    plotInternal(axs[1,0], t, Te, ylabel=r'$T_{\rm e}$ (eV)', color='g', showlabel=False, label='STREAM')
    plotInternal(axs[1,1], t[1:], Lf, ylabel=r'$L_f$ (m)', color='g', showlabel=False, label='STREAM')
    plotInternal(axs[2,0], t, Ti, ylabel=r'$T_{\rm i}$ (eV)', color='g', showlabel=False, label='STREAM')
    plotInternal(axs[2,1], t[1:], tau, ylabel=r'$\tau_{\rm D}$ (s)', color='g', showlabel=False, label='STREAM')

    Ip_mat = np.genfromtxt('Data/PlasmaCurrent_DYON.csv', delimiter=',')
    t_Ip_d = Ip_mat[:, 0]
    Ip_d = Ip_mat[:, 1]
    plotInternal(axs[0, 0], t_Ip_d, Ip_d/1e3, ylabel=r'$I_{\rm p}$ (kA)', color='k', showlabel=showlabel, label='DYON')

    ne_mat = np.genfromtxt('Data/ElectronDensity_DYON.csv', delimiter=',')
    t_ne_d = ne_mat[:, 0]
    ne_d = ne_mat[:, 1]
    plotInternal(axs[0, 1], t_ne_d, ne_d/1e17, ylabel=r'$n_{\rm e}$ ($17\cdot 10^17$m$^{-3}$)', color='k', showlabel=False, label='DYON')

    Te_mat = np.genfromtxt('Data/ElectronTemperature_DYON.csv', delimiter=',')
    t_Te_d = Te_mat[:, 0]
    Te_d = Te_mat[:, 1]
    plotInternal(axs[1, 0], t_Te_d, Te_d, ylabel=r'$T_{\rm e}$ (eV)', color='k', showlabel=False, label='DYON')

    Lf_mat = np.genfromtxt('Data/ConnectionLength_DYON.csv', delimiter=',')
    t_Lf_d = Lf_mat[:, 0]
    Lf_d = Lf_mat[:, 1]
    plotInternal(axs[1,1], t_Lf_d, Lf_d, ylabel=r'$L_f$ (m)', color='k', showlabel=False, label='DYON')

    Ti_mat = np.genfromtxt('Data/IonTemperature_DYON.csv', delimiter=',')
    t_Ti_d = Ti_mat[:, 0]
    Ti_d = Ti_mat[:, 1]
    plotInternal(axs[2, 0], t_Ti_d, Ti_d, ylabel=r'$T_{\rm i}$ (eV)', color='k', showlabel=False, label='DYON')

    tau_mat = np.genfromtxt('Data/ConfinementTime_DYON.csv', delimiter=',')
    t_tau_d = tau_mat[:, 0]
    tau_d = tau_mat[:, 1]
    plotInternal(axs[2, 1], t_tau_d, tau_d, ylabel=r'$\tau_{\rm D}$ (s)', color='k', showlabel=False, label='DYON')

    for i in range(axs.shape[0]):
        for j in range(axs.shape[1]):
            axs[i,j].set_xlim([0, 0.3])
            axs[i,j].grid(True)

    axs[0,0].set_ylim([0, 230])
    axs[0,1].set_ylim([0, 15])
    axs[1,0].set_ylim([0, 80])
    axs[1,1].set_ylim([0, 1.5e4])
    axs[2,0].set_ylim([0, 80])
    axs[2,1].set_ylim([0, 0.5])

    axs[0,0].set_yticks([0, 50, 100, 150, 200])
    axs[0,1].set_yticks([0, 5, 10, 15])
    axs[1,0].set_yticks([0, 20, 40, 60, 80])
    axs[2,0].set_yticks([0, 20, 40, 60, 80])
    axs[2,1].set_yticks([0, 0.1, 0.2, 0.3, 0.4])

    axs[0,0].legend(loc='best', frameon=False, prop={'size': 10})


def drawplot2(axs, so, toffset=0):
    t = so.grid.t[:] + toffset
    showlabel = (toffset==0)

    Vp = so.other.stream.V_p[:,0] / 1e6
    Poh = -so.other.fluid.Tcold_ohmic[:,0] * Vp
    Pequi = so.other.fluid.Tcold_ion_coll[:,0] * Vp
    #Pconve = so.other.scalar.energyloss_T_cold[:,0] / 1e6
    Pconve = so.other.stream.Tcold_transport[:,0] * Vp
    Prad = so.other.fluid.Tcold_radiation[:,0] * Vp

    PequiI = so.other.stream.Wi_e_coll[:,0] * Vp
    Pcx = -so.other.stream.Wi_chargeexchange[:,0] * Vp
    Pconvi = -so.other.stream.Wi_iontransport[:,0] * Vp

    ylbl = 'MW'
    plotInternal(axs[0], t[1:], Poh, ylabel=ylbl, color='r', linestyle='--', label='Ohmic', showlabel=showlabel)
    plotInternal(axs[0], t[1:], Pequi, ylabel=ylbl, color='b', linestyle='--', label='Equilibration', showlabel=showlabel)
    plotInternal(axs[0], t[1:], Pconve, ylabel=ylbl, color='k', linestyle='--', label='Transport', showlabel=showlabel)
    plotInternal(axs[0], t[1:], Prad, ylabel=ylbl, color='m', linestyle='--', label='Radiation + ionization', showlabel=showlabel)

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
    drawplot1(axs1, so2, toffset=so1.grid.t[-1], showlabel=True)

    fig2, axs2 = plt.subplots(1, 2, figsize=(10, 4))

    drawplot2(axs2, so1)
    drawplot2(axs2, so2, toffset=so1.grid.t[-1])

    fig1.tight_layout()
    #fig1.savefig('ITERcase.pdf')
    fig2.tight_layout()
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
        prefill = 2 * 0.8e-3 / 133.32   # Pa -> Torr
        ss1 = generate(prefill=prefill, nt=40000)
        ss1.save(f'settings1{ext}.h5')
        so1 = runiface(ss1, f'output1{ext}.h5', quiet=False)

        ss2 = STREAMSettings(ss1)
        ss2.fromOutput(f'output1{ext}.h5')
        ss2.timestep.setTmax(0.3 - ss1.timestep.tmax)
        ss2.timestep.setNumberOfSaveSteps(0)
        ss2.timestep.setNt(3000)
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


