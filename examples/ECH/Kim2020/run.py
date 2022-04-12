#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
from scipy.interpolate import interp1d
import sys

sys.path.append('../../../py')

# from run import makeplots
import PlasmaParameters as Formulas
import DREAM.Settings.Equations.ColdElectronTemperature as Tcold
import DREAM.Settings.Solver as Solver

import DREAM.Settings.Atomics as Atomics
from STREAM import STREAMOutput, STREAMSettings, runiface
import STREAM.Settings.Equations.ElectricField as ElectricField
import STREAM.Settings.Equations.IonSpecies as Ions


def generate(prefill=0.27e-3, gamma=2e-3, fractionO = 0.02, fractionC = 0, Ip=2.4e3, tmax=2e-5, nt=2000):
    """
    Generate a STREAMSettings object for a simulation with the specified
    parameters.

    :param prefill: Prefill gas pressure [Pa]
    :param gamma:   Degree of ionization
    :param Vloop:   External loop voltage [V]
    :param Vloop_t: Time vector corresponding to given loop voltage [V]
    :param j0:      Initial plasma current density [A/m^2]
    :param tmax:    Simulation time [s]
    :param nt:      Number of time steps
    """
    n0 = 4.8e20 * prefill  # Initial total deuterium density
    nD = n0 * np.array([[1 - gamma], [gamma]])
    nO = fractionO * n0
    if fractionC <= 0:
        nC = 1e3
    else:
        nC = fractionC * n0

    a = 0.5  # Plasma minor radius [m]
    R0 = 1.8  # Plasma major radius [m]
    Btor = 2.7  # * 6.2 / R0  # Toroidal magnetic field [T]
    V_vessel = 56  # Vacuum vessel volume

    Te0 = 1  # electron temperature [eV]
    Ti0 = 0.026  # ion temperature [eV]

    # Initial electric field
    j0 = Ip / (a ** 2 * np.pi)
    E0 = j0 / Formulas.evaluateSpitzerConductivity(n=nD[1], T=Te0, Z=1)

    P_inj = 800e3
    f_o = 0
    f_x = 1.0
    theta = 14 * np.pi / 180
    phi = 20 * np.pi / 180
    N = 2

    t_start = 0
    t_dur = 0.26
    t_end = t_start + t_dur
    Nt = 53
    t = np.linspace(t_start, t_end, Nt)

    c1 = 1.55
    c2 = 0.52
    c3 = 0.1
    Y_DD = c1 - c2 * (1 - np.exp(-(t - t_start) / c3))

    t_DC_vec = np.array([0.04, 0.205, 0.21, 0.215, 0.22, 0.225, 0.30]) - 0.04
    Y_DC_vec = np.array([0.50, 0.500, 0.51, 0.535, 0.58, 0.645, 2.50])

    Y_DC_fun = interp1d(t_DC_vec, Y_DC_vec, kind='linear')
    Y_DC = Y_DC_fun(t) / 100

    t_V_vec1 = np.array([0.04, 0.05, 0.08, 0.10, 0.12]) - 0.04
    V_vec1 = np.array([3.6, 3.90, 4.35, 4.45, 4.49])
    t_V_vec2 = np.array([0.12, 0.127, 0.15, 0.20, 0.25, 0.31]) - 0.04
    V_vec2 = np.array([4.49, 4.000, 3.10, 2.44, 2.31, 2.30])

    V_fun1 = interp1d(t_V_vec1, V_vec1, kind='cubic')
    V_fun2 = interp1d(t_V_vec2, V_vec2, kind='cubic')
    V_loop = np.append(V_fun1(t[:16]), V_fun2(t[16:]))

    '''
    t_plot = np.linspace(0.04, 0.3, Nt)
    plt.plot(t_plot, Y_DD, 'k')
    plt.plot(t_plot, Y_DC * 100, 'r')
    plt.xlim([0, 0.3])
    plt.ylim([0, 3])
    plt.grid(True)
    plt.show()
    plt.plot(t_plot, V_loop, 'k')
    plt.xlim([0, 0.3])
    plt.ylim([0, 5])
    plt.grid(True)
    plt.show()
    #'''

    ss = STREAMSettings()

    ss.atomic.adas_interpolation = Atomics.ADAS_INTERP_BILINEAR

    # Electric field
    ss.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
    ss.eqsys.E_field.setInitialProfile(E0)
    Lp = float(scipy.constants.mu_0 * R0 * (np.log(8 * R0 / a) + 0.25 - 2))
    ss.eqsys.E_field.setInductances(Lp=Lp, Lwall=9.1e-6, M=2.49e-6, Rwall=1e6)
    ss.eqsys.E_field.setCircuitVloop(V_loop, t)

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

    # Recycling coefficients
    iD = ss.eqsys.n_i.getIndex('D')
    ss.eqsys.n_i.ions[iD].setRecyclingCoefficient('D', Y_DD, t)
    ss.eqsys.n_i.ions[iD].setRecyclingCoefficient('C', Y_DC, t)
    
    iC = ss.eqsys.n_i.getIndex('C')
    ss.eqsys.n_i.ions[iC].setRecyclingCoefficient('C', 0.015)

    iO = ss.eqsys.n_i.getIndex('O')
    ss.eqsys.n_i.ions[iO].setRecyclingCoefficient('C', 1)
    ss.eqsys.n_i.ions[iO].setRecyclingCoefficient('O', 1)

    # Radial grid
    ss.radialgrid.setB0(Btor)
    ss.radialgrid.setMinorRadius(a)
    ss.radialgrid.setMajorRadius(R0)
    ss.radialgrid.setWallRadius(a)
    ss.radialgrid.setVesselVolume(V_vessel)
    ss.radialgrid.setIref(10e3)
    ss.radialgrid.setBv(2.2e-3)
    ss.radialgrid.setConnectionLengthFactor(1.0)

    ss.radialgrid.setECHParameters(P_inj, f_o, f_x, theta, phi, N)

    # Disable kinetic grids
    ss.hottailgrid.setEnabled(False)
    ss.runawaygrid.setEnabled(False)

    # Numerical settings
    ss.solver.setType(Solver.NONLINEAR)
    ss.solver.preconditioner.setEnabled(False)
    ss.timestep.setTmax(tmax)
    ss.timestep.setNt(nt)
    ss.timestep.setNumberOfSaveSteps(10000)

    ss.other.include('fluid', 'stream', 'scalar', 'fluid/Tcold_ECH')

    return ss


def drawplot1(axs, so, toffset=0, showlabel=False, save=False, first=True):
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
        if first:
            worab = 'w'
        else:
            worab = 'ab'
        t_csv = open('Data/time_STREAM.csv', worab)
        np.savetxt(t_csv, t)
        t_csv.close()
        Ip_csv = open('Data/PlasmaCurrent_STREAM.csv', worab)
        np.savetxt(Ip_csv, Ip)
        Ip_csv.close()
        ne_csv = open('Data/ElectronDensity_STREAM.csv', worab)
        np.savetxt(ne_csv, ne)
        ne_csv.close()
        Te_csv = open('Data/ElectronTemperature_STREAM.csv', worab)
        np.savetxt(Te_csv, Te)
        Te_csv.close()
        Ti_csv = open('Data/IonTemperature_STREAM.csv', worab)
        np.savetxt(Ti_csv, Ti)
        Ti_csv.close()
        Lf_csv = open('Data/ConnectionLength_STREAM.csv', worab)
        np.savetxt(Lf_csv, Lf)
        Lf_csv.close()
        tau_csv = open('Data/ConfinementTime_STREAM.csv', worab)
        np.savetxt(tau_csv, tau)
        tau_csv.close()

    plotInternal(axs[0, 0], t, Ip / 1e6, ylabel=r'$I_{\rm p}$ (kA)', color='g', showlabel=showlabel, label='STREAM')
    plotInternal(axs[0, 1], t, ne / 1e19, ylabel=r'$n_{\rm e}$ ($1\cdot 10^{19}$m$^{-3}$)', color='g', showlabel=False,
                 label='STREAM')
    plotInternal(axs[1, 0], t, Te, ylabel=r'$T_{\rm e}$ (eV)', color='g', showlabel=False, label='STREAM')
    plotInternal(axs[1, 1], t[1:], Lf, ylabel=r'$L_f$ (m)', color='g', showlabel=False, label='STREAM', log=True)
    plotInternal(axs[2, 0], t, Ti, ylabel=r'$T_{\rm i}$ (eV)', color='g', showlabel=False, label='STREAM')
    plotInternal(axs[2, 1], t[1:], tau, ylabel=r'$\tau_{\rm D}$ (s)', color='g', showlabel=False, label='STREAM')

    for i in range(axs.shape[0]):
        for j in range(axs.shape[1]):
            axs[i, j].set_xlim([0, 0.3])
            axs[i, j].grid(True)

    axs[0, 0].set_ylim([0, 0.3])
    axs[0, 1].set_ylim([0, 1.1])
    axs[1, 0].set_ylim([0, 350])
    axs[1, 1].set_ylim([1e0, 1e11])
    axs[2, 0].set_ylim([0, 100])
    axs[2, 1].set_ylim([0, 0.04])
    '''
    axs[0, 0].set_yticks([0, 50, 100, 150, 200])
    axs[0, 1].set_yticks([0, 5, 10, 15])
    axs[1, 0].set_yticks([0, 20, 40, 60, 80])
    axs[2, 0].set_yticks([0, 20, 40, 60, 80])
    axs[2, 1].set_yticks([0, 0.1, 0.2, 0.3, 0.4])

    axs[0, 0].legend(loc='best', frameon=False, prop={'size': 10})
    '''


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


def makeplots(so1, so2, toffset):
    fig1, axs1 = plt.subplots(3, 2, figsize=(7, 10))

    drawplot1(axs1, so1, toffset=toffset)
    drawplot1(axs1, so2, toffset=toffset + so1.grid.t[-1], showlabel=True, first=False)


    fig1.tight_layout()
    plt.show()


def main(argv):
    FONTSIZE = 16
    plt.rcParams.update({'font.size': FONTSIZE})

    parser = argparse.ArgumentParser()

    parser.add_argument('-e', '--extension', help="Append the specified string to the end of file names",
                        dest="extension", action='store', default='')
    parser.add_argument('-n', '--no-plot', help="Only run simulations, don't generate plots", dest="plot",
                        action='store_false', default=True)
    parser.add_argument('-s', '--skip', help="Skip the simulation and load output from already existing files",
                        action='store_true', dest="skip", default=False)

    settings = parser.parse_args()

    ext = '' if not settings.extension else '_' + settings.extension

    if not settings.skip:
        ss1 = generate(nt=100000)
        ss1.save(f'settings1{ext}.h5')
        so1 = runiface(ss1, f'output1{ext}.h5', quiet=False)

        ss2 = STREAMSettings(ss1)
        ss2.fromOutput(f'output1{ext}.h5')
        ss2.timestep.setTmax(0.3 - ss1.timestep.tmax)
        ss2.timestep.setNumberOfSaveSteps(0)
        ss2.timestep.setNt(20000)
        ss2.save(f'settings2{ext}.h5')
        so2 = runiface(ss2, f'output2{ext}.h5', quiet=False)
    else:
        so1 = STREAMOutput(f'output1{ext}.h5')
        so2 = STREAMOutput(f'output2{ext}.h5')

    if settings.plot:
        makeplots(so1, so2, toffset=0.04)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


