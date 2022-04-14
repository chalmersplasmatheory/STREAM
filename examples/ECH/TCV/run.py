#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
from scipy.interpolate import interp1d
import sys
import h5py
from scipy.signal import savgol_filter

sys.path.append('../../../py')

# from run import makeplots
import PlasmaParameters as Formulas
import DREAM.Settings.Equations.ColdElectronTemperature as Tcold
import DREAM.Settings.Solver as Solver
import DREAM.Settings.Equations.RunawayElectrons as Runaways

import DREAM.Settings.Atomics as Atomics
from STREAM import STREAMOutput, STREAMSettings, runiface
import STREAM.Settings.Equations.ElectricField as ElectricField
import STREAM.Settings.Equations.IonSpecies as Ions


def generate(gamma=2e-3, fractionNe = 0.02, tmax=1e-5, nt=2000, tstart=-0.01, tend=0.11):
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
    n0 = 5e17  # Initial total deuterium density # ?
    nD = n0 * np.array([[1 - gamma], [gamma]])
    if fractionNe <= 0: # fractionNe ?
        nNe = 1e3
    else:
        nNe = fractionNe * n0

    hf = h5py.File('TCV65108.h5', 'r')
    t_loop = np.array(hf.get('Loop voltage').get('x'))
    V_loop_osc =- np.array(hf.get('Loop voltage').get('z'))
    V_loop = savgol_filter(V_loop_osc, 105, 3)
    t_fluxD = np.array(hf.get('Particle flux (D2)').get('x'))
    fluxD = np.array(hf.get('Particle flux (D2)').get('z'))
    t_Ip = np.array(hf.get('Plasma current').get('x'))
    Ip = np.array(hf.get('Plasma current').get('z'))
    B = np.array(hf.get('Toroidal magnetic field').get('z'))


    R0 = 0.89  # Plasma major radius [m]
    Btor = 1.45   # Toroidal magnetic field [T]
    V_vessel = 4.632  # Vacuum vessel volume

    t_Vp_osc = np.array(hf.get('Plasma volume').get('x'))
    V_p_osc = np.array(hf.get('Plasma volume').get('z'))
    V_p_s = savgol_filter(V_p_osc, 75, 3)
    V_initfun = interp1d(np.append(np.array([-0.02]), t_Vp_osc), np.append(np.array([V_vessel]), V_p_osc), 'cubic')
    t_Vp = np.linspace(-0.02, t_Vp_osc[-1])
    V_p = V_initfun(t_Vp)
    a = np.sqrt(V_p / (2 * np.pi ** 2 * R0))

    hf.close()

    Te0 = 1  # electron temperature [eV]
    Ti0 = 0.026  # ion temperature [eV]

    # Initial electric field
    it0 = int((tstart - t_Ip[0])/(t_Ip[-1] - t_Ip[0]) * t_Ip.shape[0])
    itend = int((tend - t_Ip[0]) / (t_Ip[-1] - t_Ip[0]) * t_Ip.shape[0])
    it0a = np.where(np.min(np.abs(t_Vp - tstart)))
    j0 = Ip[it0] / (a[it0a] ** 2 * np.pi)
    E0 = j0 / Formulas.evaluateSpitzerConductivity(n=nD[1], T=Te0, Z=1)

    P_inj = 675e3 # 650e3 - 700e3
    f_o = 0 # ??
    f_x = 1.0 # ??
    theta = 10 * np.pi / 180 # ??
    phi = 50 * np.pi / 180 # ??
    N = 2 # ??

    # Impurities ??

    ss = STREAMSettings()

    ss.atomic.adas_interpolation = Atomics.ADAS_INTERP_BILINEAR

    # Electric field
    ss.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
    ss.eqsys.E_field.setInitialProfile(E0)
    Lp = float(scipy.constants.mu_0 * R0 * (np.log(8 * R0 / 0.25) + 0.25 - 2))
    ss.eqsys.E_field.setInductances(Lp=Lp, Lwall=9.1e-6, M=2.49e-6, Rwall=1e6) # ?
    ss.eqsys.E_field.setCircuitVloop(V_loop, t_loop)

    # Electron temperature
    ss.eqsys.T_cold.setType(Tcold.TYPE_SELFCONSISTENT)
    ss.eqsys.T_cold.setInitialProfile(Te0)
    ss.eqsys.T_cold.setRecombinationRadiation(recombination=Tcold.RECOMBINATION_RADIATION_INCLUDED)

    # Ions
    ss.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=nD, r=np.array([0]), T=Ti0)
    #ss.eqsys.n_i.addIon(name='Ne', Z=10, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=nNe, r=np.array([0]), T=Ti0)

    # Disable runaway
    ss.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
    ss.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)

    # Recycling coefficients
    t = np.linspace(tstart, tend)
    c1 = 1.022
    c2 = 0.02
    c3 = 0.1
    Y_DD = c1 - c2 * (1 - np.exp(-(t-tstart) / c3))
    iD = ss.eqsys.n_i.getIndex('D')
    ss.eqsys.n_i.ions[iD].setRecyclingCoefficient('D', 1) # ?

    it0 = int((tstart - t_fluxD[0]) / (t_fluxD[-1] - t_fluxD[0]) * t_fluxD.shape[0])
    itend = int((tend - t_fluxD[0]) / (t_fluxD[-1] - t_fluxD[0]) * t_fluxD.shape[0])

    plt.show()
    ss.eqsys.n_i.setFueling('D', fluxD[it0:itend]*1.5e-1, times=t_fluxD[it0:itend]) # ?

    # Radial grid
    ss.radialgrid.setB0(Btor)
    ss.radialgrid.setMinorRadius(a, t_Vp)
    ss.radialgrid.setMajorRadius(R0)
    ss.radialgrid.setWallRadius(0.25)
    ss.radialgrid.setVesselVolume(V_vessel)
    ss.radialgrid.setIref(5e3) # Maximal plasma current 1.2 MA, KSTAR 2 MA so same I_ref? Half? B is half for TCV ??
    ss.radialgrid.setBv(2.0e-3) # ??
    ss.radialgrid.setConnectionLengthFactor(1.0) # ??

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

    hf = h5py.File('TCV65108.h5', 'r')
    t_Ip_d = np.array(hf.get('Plasma current').get('x')) + 0.02
    Ip_d = np.array(hf.get('Plasma current').get('z'))
    hf.close()
    hf = h5py.File('firdens.h5', 'r')
    t_FIR = np.array(hf.get('FIR').get('x')) + 0.02
    FIR = np.array(hf.get('FIR').get('z'))
    hf.close()

    plotInternal(axs[0], t_Ip_d, Ip_d / 1e6, ylabel=r'$I_{\rm p}$ (MA)', color='tab:red', showlabel=showlabel, label='STREAM')
    plotInternal(axs[0], t, Ip / 1e6, ylabel=r'$I_{\rm p}$ (MA)', color='tab:blue', showlabel=showlabel, label='STREAM')
    plotInternal(axs[1], t_FIR, FIR / 1e19, ylabel=r'$n_{\rm e}$ ($1\cdot 10^{19}$m$^{-3}$)', color='tab:red', showlabel=False,
                 label='STREAM')
    plotInternal(axs[1], t, ne / 1e19, ylabel=r'$n_{\rm e}$ ($1\cdot 10^{19}$m$^{-3}$)', color='tab:blue', showlabel=False,
                 label='STREAM')
    plotInternal(axs[2], t, Te, ylabel=r'$T_{\rm e}$ (eV)', color='tab:blue', showlabel=False, label='STREAM')
    '''
    plotInternal(axs[1, 1], t[1:], Lf, ylabel=r'$L_f$ (m)', color='tab:blue', showlabel=False, label='STREAM', log=True)
    plotInternal(axs[2, 0], t, Ti, ylabel=r'$T_{\rm i}$ (eV)', color='tab:blue', showlabel=False, label='STREAM')
    plotInternal(axs[2, 1], t[1:], tau, ylabel=r'$\tau_{\rm D}$ (s)', color='tab:blue', showlabel=False, label='STREAM')
    
    for i in range(axs.shape[0]):
        for j in range(axs.shape[1]):
            #axs[i, j].set_xlim([0, 0.3])
            #axs[i, j].grid(True)
    '''
    axs[0].set_ylim([0, 0.4])
    axs[1].set_ylim([0, 2])
    #axs[2, 0].set_ylim([0, 100])


def drawplot2(axs, so, toffset=0, showlabel=False, save=False, first=True):
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

    plotInternal(axs[0, 0], t, Ip / 1e3, ylabel=r'$I_{\rm p}$ (kA)', color='g', showlabel=showlabel, label='STREAM')
    plotInternal(axs[0, 1], t, ne / 1e19, ylabel=r'$n_{\rm e}$ ($1\cdot 10^{19}$m$^{-3}$)', color='g', showlabel=False,
                 label='STREAM')
    plotInternal(axs[1, 0], t, Te, ylabel=r'$T_{\rm e}$ (eV)', color='g', showlabel=False, label='STREAM')
    plotInternal(axs[1, 1], t[1:], Lf, ylabel=r'$L_f$ (m)', color='g', showlabel=False, label='STREAM', log=True)
    plotInternal(axs[2, 0], t, Ti, ylabel=r'$T_{\rm i}$ (eV)', color='g', showlabel=False, label='STREAM')
    plotInternal(axs[2, 1], t[1:], tau, ylabel=r'$\tau_{\rm D}$ (s)', color='g', showlabel=False, label='STREAM')

    for i in range(axs.shape[0]):
        for j in range(axs.shape[1]):
            #axs[i, j].set_xlim([0, 0.3])
            axs[i, j].grid(True)

    #axs[0, 0].set_ylim([0, 0.3])
    #axs[0, 1].set_ylim([0, 1.1])
    #axs[1, 0].set_ylim([0, 350])
    #axs[1, 1].set_ylim([1e0, 1e11])
    #axs[2, 0].set_ylim([0, 100])
    #axs[2, 1].set_ylim([0, 0.04])
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
    fig1, axs1 = plt.subplots(3, 1, figsize=(7, 10))

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
        ss1 = generate(nt=10000)
        ss1.save(f'settings1{ext}.h5')
        so1 = runiface(ss1, f'output1{ext}.h5', quiet=False)

        ss2 = STREAMSettings(ss1)
        ss2.fromOutput(f'output1{ext}.h5')
        ss2.timestep.setTmax(0.12 - ss1.timestep.tmax)
        ss2.timestep.setNumberOfSaveSteps(0)
        ss2.timestep.setNt(50000)
        ss2.save(f'settings2{ext}.h5')
        so2 = runiface(ss2, f'output2{ext}.h5', quiet=False)
    else:
        so1 = STREAMOutput(f'output1{ext}.h5')
        so2 = STREAMOutput(f'output2{ext}.h5')

    if settings.plot:
        makeplots(so1, so2, toffset=0)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


