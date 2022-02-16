#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
from scipy.interpolate import interp1d
import sys

sys.path.append('../../../py')

#from run import makeplots
import PlasmaParameters as Formulas
import DREAM.Settings.Equations.ColdElectronTemperature as Tcold
import DREAM.Settings.Solver as Solver
import DREAM.Settings.Equations.ElectricField as ElectricField_D
import DREAM.Settings.Atomics as Atomics
from STREAM import STREAMOutput, STREAMSettings, runiface
import STREAM.Settings.Equations.ElectricField as ElectricField
import STREAM.Settings.Equations.IonSpecies as Ions


def generate(prefill=5e-5, gamma=2e-3, fractionO = 0.001, fractionC = 0, Ip=2.4e3, tmax=1e-4, nt=2000, selfconsistent = False):
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
    #n0 = 1.296e18
    nD = n0 * np.array([[1-gamma], [gamma]])
    nO = fractionO * n0
    nC = fractionC * n0

    Btor = 2.4      # Toroidal magnetic field [T]
    R0 = 2.96       # Plasma major radius [m]
    r_wall = 1      # Distance between plasma centre and passive structure [m]
    #r_wall = 1.3
    V_vessel = 100  # Vacuum vessel volume

    t = np.linspace(0, 0.3, 100)
    t_a = np.array([0, 0.017, 0.05, 0.085, 0.14, 0.19, 0.25, 0.3])
    V_p = np.array([100, 80, 56, 48, 52, 51.75, 54.25, 56])
    a_vec = np.sqrt(V_p / (2 * np.pi ** 2 * R0))
    a_fun = interp1d(t_a, a_vec, kind='cubic')
    a = a_fun(t)  # Plasma minor radius [m]

    tVloop = [0, 0.02, 0.0325, 0.0475, 0.08, 0.1, 0.125, 0.13, 0.15, 0.20, 0.22, 0.23, 0.25, 0.3, 0.335, 0.35, 0.37,
               0.4, 0.45, 0.5]
    Vloop = [11, 21.25, 26, 26.25, 24, 16.5, 8.25, 7.9, 7.75, 7.5, 7.25, 6.5, 6.5, 6.75, 6.75, 6, 4.75, 4.25, 4.5,
               3.60]

    c1 = 1.1
    #c2 = 0.09
    c2 = 0.05
    c3 = 0.1

    Te0 = 1     # electron temperature [eV]
    Ti0 = 0.026 # ion temperature [eV]

    # Initial electric field
    sigma = Formulas.evaluateSpitzerConductivity(n=nD[1], T=Te0, Z=1)
    j0 = Ip / (a_vec[0] ** 2 * np.pi)  # 943?
    E0 = j0 / sigma

    ss = STREAMSettings()

    ss.atomic.adas_interpolation = Atomics.ADAS_INTERP_BILINEAR

    # Electric field
    if selfconsistent:
        R = 7.5e-4  # Resistance in MK2 structure [Ohm]
        L = 9.1e-6  # Inductance in MK2 structure [H]
        wall_time = L / R
        ss.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
        ss.eqsys.E_field.setInitialProfile(efield=E0)
        ss.eqsys.E_field.setBoundaryCondition(ElectricField_D.BC_TYPE_TRANSFORMER,
                                                       V_loop_wall_R0=Vloop / R0, times=t,
                                                       inverse_wall_time=1 / wall_time, R0=R0)
    else:
        ss.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
        ss.eqsys.E_field.setInitialProfile(efield=E0)
        ss.eqsys.E_field.setInductances(Lp=5.19e-6, Lwall=9.1e-6, M=2.49e-6, Rwall=7.5e-4)
        ss.eqsys.E_field.setCircuitVloop(Vloop, tVloop)

    # Electron temperature
    ss.eqsys.T_cold.setType(Tcold.TYPE_SELFCONSISTENT)
    ss.eqsys.T_cold.setInitialProfile(Te0)
    ss.eqsys.T_cold.setRecombinationRadiation(recombination=Tcold.RECOMBINATION_RADIATION_INCLUDED)

    # Ions
    ss.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=nD, r=np.array([0]), T=Ti0)
    ss.eqsys.n_i.addIon(name='C', Z=6, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=nC, r=np.array([0]), T=Ti0)
    ss.eqsys.n_i.addIon(name='O', Z=8, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=nO, r=np.array([0]), T=Ti0)

    # Disable runaway
    ss.eqsys.n_re.setAvalanche(False)
    ss.eqsys.n_re.setDreicer(False)

    # Recycling coefficients (unused)
    #ss.eqsys.n_i.setJET_CWrecycling()
    for ion in ss.eqsys.n_i.ions:
        if ion.name == 'D':
            ion.setRecyclingCoefficient('C', 0.015)
        elif ion.name == 'O':
            ion.setRecyclingCoefficient('C', 1.0)
            ion.setRecyclingCoefficient('O', 1.0)

    # Radial grid
    ss.radialgrid.setB0(Btor)
    ss.radialgrid.setMinorRadius(a, t)
    ss.radialgrid.setMajorRadius(R0)
    ss.radialgrid.setWallRadius(r_wall)
    ss.radialgrid.setVesselVolume(V_vessel)

    ss.radialgrid.setRecyclingCoefficient1(c1)
    ss.radialgrid.setRecyclingCoefficient2(c2)
    ss.radialgrid.setRecyclingCoefficient3(c3)
    
    # Disable kinetic grids
    ss.hottailgrid.setEnabled(False)
    ss.runawaygrid.setEnabled(False)

    # Numerical settings
    ss.solver.setType(Solver.NONLINEAR)
    ss.solver.preconditioner.setEnabled(False)
    ss.timestep.setTmax(tmax)
    ss.timestep.setNt(nt)

    ss.other.include('fluid', 'stream', 'scalar')

    return ss


def drawplot1(axs, so, toffset=0, showlabel=False, save=False):
    """
    Draw a plot with a output from the given STREAMOutput object.
    """
    t = so.grid.t[:] + toffset
    Ip = so.eqsys.I_p[:, 0]
    ne = so.eqsys.n_cold[:, 0]
    Te = so.eqsys.T_cold[:, 0]
    Ti = so.eqsys.W_i.getTemperature()['D'][:, 0]
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

    plotInternal(axs[0,0], t, Ip, ylabel=r'$I_{\rm p}$ (A)', color='g', showlabel=showlabel, label='STREAM')
    plotInternal(axs[0,1], t, ne, ylabel=r'$n_{\rm e}$ (m$^{-3}$)', color='g', showlabel=showlabel, label='STREAM')
    plotInternal(axs[1,0], t, Te, ylabel=r'$T_{\rm e}$ (eV)', color='g', showlabel=showlabel, label='STREAM')
    plotInternal(axs[1,1], t[1:], Lf, ylabel=r'$L_f$ (m)', color='g', showlabel=False, label='STREAM', log=True)
    plotInternal(axs[2,0], t, Ti, ylabel=r'$T_{\rm i}$ (eV)', color='g', showlabel=showlabel, label='STREAM')
    plotInternal(axs[2,1], t[1:], tau, ylabel=r'$\tau_{\rm D}$ (s)', color='g', showlabel=showlabel, label='STREAM')

    Ip_mat = np.genfromtxt('Data/PlasmaCurrent_DYON.csv', delimiter=',')
    t_Ip_d = Ip_mat[:,0]
    Ip_d = Ip_mat[:, 1]
    plotInternal(axs[0, 0], t_Ip_d, Ip_d, ylabel=r'$I_{\rm p}$ (A)', color='k', showlabel=showlabel, label='DYON')

    Ip_mat_m = np.genfromtxt('Data/PlasmaCurrent_measured.csv', delimiter=',')
    t_Ip_m = Ip_mat_m[:, 0]
    Ip_m = Ip_mat_m[:, 1]
    plotInternal(axs[0, 0], t_Ip_m, Ip_m, ylabel=r'$I_{\rm p}$ (A)', color='silver', showlabel=showlabel, label='Measured')

    ne_mat = np.genfromtxt('Data/ElectronDensity_DYON.csv', delimiter=',')
    t_ne_d = ne_mat[:, 0]
    ne_d = ne_mat[:, 1]
    plotInternal(axs[0, 1], t_ne_d, ne_d, ylabel=r'$n_{\rm e}$ (m$^{-3}$)', color='k', showlabel=showlabel, label='DYON')

    Te_mat = np.genfromtxt('Data/ElectronTemperature_DYON.csv', delimiter=',')
    t_Te_d = Te_mat[:, 0]
    Te_d = Te_mat[:, 1]
    plotInternal(axs[1, 0], t_Te_d, Te_d, ylabel=r'$T_{\rm e}$ (eV)', color='k', showlabel=showlabel, label='DYON')

    Lf_mat = np.genfromtxt('Data/ConnectionLength_DYON.csv', delimiter=',')
    t_Lf_d = Lf_mat[:, 0]
    Lf_d = Lf_mat[:, 1]
    plotInternal(axs[1,1], t_Lf_d, Lf_d, ylabel=r'$L_f$ (m)', color='k', showlabel=False, label='DYON', log=True)

    Ti_mat = np.genfromtxt('Data/IonTemperature_DYON.csv', delimiter=',')
    t_Ti_d = Ti_mat[:, 0]
    Ti_d = Ti_mat[:, 1]
    plotInternal(axs[2, 0], t_Ti_d, Ti_d, ylabel=r'$T_{\rm i}$ (eV)', color='k', showlabel=showlabel, label='DYON')

    tau_mat = np.genfromtxt('Data/ConfinementTime_DYON.csv', delimiter=',')
    t_tau_d = tau_mat[:, 0]
    tau_d = tau_mat[:, 1]
    plotInternal(axs[2, 1], t_tau_d, tau_d, ylabel=r'$\tau_{\rm D}$ (s)', color='k', showlabel=showlabel, label='DYON')


    for i in range(axs.shape[0]):
        for j in range(axs.shape[1]):
            axs[i,j].set_xlim([0, 0.3])
            axs[i,j].grid(True)

    axs[0,0].set_ylim([0, 6e5])
    axs[0,1].set_ylim([0, 6e18])
    axs[1,0].set_ylim([0, 400])
    #axs[1, 1].set_ylim([1, 1.5e4])
    axs[2,0].set_ylim([0, 400])
    axs[2,1].set_ylim([0, 0.1])

    axs[0,0].set_yticks([0, 2e5, 4e5, 6e5])
    axs[0,1].set_yticks([0, 2e18, 4e18, 6e18])
    axs[1,0].set_yticks([0, 100, 200, 300, 400])
    axs[2,0].set_yticks([0, 100, 200, 300, 400])
    axs[2,1].set_yticks([0, 0.05, 0.1])

    axs[0,0].legend(loc='best', frameon=False, prop={'size': 10})

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


    fig1.tight_layout()
    #fig1.savefig('JETcase.pdf')
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
        prefill = 2 * 2.7e-3 / 133.32   # Pa -> Torr
        ss1 = generate(prefill=prefill)
        ss1.save(f'settings1{ext}.h5')
        so1 = runiface(ss1, f'output1{ext}.h5', quiet=False)

        ss2 = STREAMSettings(ss1)
        ss2.fromOutput(f'output1{ext}.h5')
        ss2.timestep.setTmax(0.3 - ss1.timestep.tmax)
        ss2.timestep.setNt(5000)
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


