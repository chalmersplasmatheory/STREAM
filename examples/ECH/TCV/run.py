#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
from scipy.interpolate import interp1d, interp2d
import sys
import h5py
from scipy.signal import savgol_filter
from scipy.optimize import least_squares
import matplotlib as mpl

FONTSIZE = 20
mpl.rcParams.update({'text.usetex': True, 'font.family': 'sans', 'font.size': FONTSIZE})

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

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def generate(gamma=2e-3, Z_eff = 3, tmax=1e-5, nt=2000, tstart=-0.01, tend=0.4, connectionLengthFactor=3.0):
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
    fractionC = (Z_eff - 1) / (36 - 6*Z_eff)
    n0 = 5e17  # Initial total deuterium density # ?
    nD = n0 * np.array([[1 - gamma], [gamma]])
    if fractionC <= 0:
        nC = 1e3
    else:
        nC = fractionC * n0
    nC = 1e3
    nNe = 1e3

    hf = h5py.File('TCV65108_old2.h5', 'r')

    t_loop = np.array(hf.get('Loop voltage').get('x'))
    V_loop_osc =- np.array(hf.get('Loop voltage').get('z'))
    it0_loop = int((tstart - t_loop[0]) / (t_loop[-1] - t_loop[0]) * t_loop.shape[0])-1
    itmin_loop = int((0.05 - t_loop[0]) / (t_loop[-1] - t_loop[0]) * t_loop.shape[0])-1
    itend_loop = int((0.5 - t_loop[0]) / (t_loop[-1] - t_loop[0]) * t_loop.shape[0])-1
    t_loop = t_loop[it0_loop:itend_loop]-tstart
    V_loop = np.append(savgol_filter(V_loop_osc[it0_loop:itmin_loop], 205, 3), savgol_filter(V_loop_osc[itmin_loop:itend_loop], 505, 3))#*1.2)
    #V_loop *= np.linspace(1, 1.2, V_loop.shape[0])
    #plt.plot(t_loop, V_loop)
    #plt.xlim([0,0.4])
    #plt.ylim([-1, 10])
    #plt.show()

    t_fluxD = np.array(hf.get('Particle flux (D2)').get('x'))
    fluxD = np.array(hf.get('Particle flux (D2)').get('z'))
    it0_flux = int((tstart - t_fluxD[0]) / (t_fluxD[-1] - t_fluxD[0]) * t_fluxD.shape[0])-1
    t_fluxD = t_fluxD[it0_flux:]-tstart
    fluxD = fluxD[it0_flux:]

    t_Ip = np.array(hf.get('Plasma current').get('x'))
    Ip = np.array(hf.get('Plasma current').get('z'))
    it0_Ip = int((tstart - t_Ip[0]) / (t_Ip[-1] - t_Ip[0]) * t_Ip.shape[0])-1
    t_Ip = t_Ip[it0_Ip:] - tstart
    Ip = Ip[it0_Ip:]

    #B = np.array(hf.get('Toroidal magnetic field').get('z'))


    R0 = 0.89  # Plasma major radius [m]
    Btor = 1.45   # Toroidal magnetic field [T]
    V_vessel = 4.632  # Vacuum vessel volume

    t_Vp_osc = np.array(hf.get('Plasma volume').get('x'))
    V_p_osc = np.array(hf.get('Plasma volume').get('z'))
    V_p_s = savgol_filter(V_p_osc, 75, 3)
    V_initfun = interp1d(np.append(np.array([-0.02-tstart]), t_Vp_osc-tstart), np.append(np.array([V_vessel]), V_p_s), 'cubic')
    t_Vp = np.linspace(0, t_Vp_osc[-1]-tstart)
    V_p = V_initfun(t_Vp)
    a = np.sqrt(V_p / (2 * np.pi ** 2 * R0))
    #plt.plot(t_Vp, a)
    #plt.show()
    hf.close()

    l_i = 0.0

    Te0 = 1  # electron temperature [eV]
    Ti0 = 0.026  # ion temperature [eV]

    # Initial electric field
    j0 = Ip[0] / (a[0] ** 2 * np.pi)
    E0 = j0 / Formulas.evaluateSpitzerConductivity(n=nD[1], T=Te0, Z=1)

    P_inj = 665e3
    f_o = 0
    f_x = 1.0
    theta = 10 * np.pi / 180
    phi = 90.0 * np.pi / 180 #np.pi/2 # ??
    N = 2

    # Impurities ??

    ss = STREAMSettings()

    ss.atomic.adas_interpolation = Atomics.ADAS_INTERP_BILINEAR

    # Electric field
    ss.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
    ss.eqsys.E_field.setInitialProfile(E0)
    Lp = float(scipy.constants.mu_0 * R0 * (np.log(8 * R0 / 0.25) + l_i/2 - 2))
    ss.eqsys.E_field.setInductances(Lp=Lp, Lwall=9.1e-6, M=2.49e-6, Rwall=1e6) # ?
    ss.eqsys.E_field.setCircuitVloop(V_loop, t_loop)

    # Electron temperature
    ss.eqsys.T_cold.setType(Tcold.TYPE_SELFCONSISTENT)
    ss.eqsys.T_cold.setInitialProfile(Te0)
    ss.eqsys.T_cold.setRecombinationRadiation(recombination=Tcold.RECOMBINATION_RADIATION_INCLUDED)

    # Ions
    ss.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=nD, r=np.array([0]), T=Ti0)
    ss.eqsys.n_i.addIon(name='C', Z=6, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=nC, r=np.array([0]), T=Ti0)
    ss.eqsys.n_i.addIon(name='Ne', Z=10, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=nNe, r=np.array([0]), T=Ti0)

    # Disable runaway
    ss.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
    ss.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)

    # Recycling coefficients
    t = np.linspace(0, tend-tstart)
    c1 = 1.022
    c2 = 0.02
    c3 = 0.1
    Y_DD = c1 - c2 * (1 - np.exp(-(t) / c3))
    iD = ss.eqsys.n_i.getIndex('D')
    ss.eqsys.n_i.ions[iD].setRecyclingCoefficient('D', 1) # ?

    iC = ss.eqsys.n_i.getIndex('C')
    ss.eqsys.n_i.ions[iD].setRecyclingCoefficient('C', [0.015, 0.05, 0.015, 0.010], [0, 0.00001, 0.05, 0.4])  # ?
    iNe = ss.eqsys.n_i.getIndex('Ne')
    ss.eqsys.n_i.ions[iNe].setRecyclingCoefficient('Ne', 1)  # ?

    g = gaussian(t_fluxD, 0.18, 0.03)
    #plt.plot(t_fluxD, g)
    #plt.xlim([0, 0.41])
    #plt.show()
    Fuel_Ne = 5e18
    diff = 1e-6
    F_Ne_fun = interp1d(np.array([-1,  -tstart-diff, -tstart, 10]), np.array([0, 0, Fuel_Ne, Fuel_Ne]))
    F_Ne = F_Ne_fun(t_fluxD)
    #plt.plot(np.array([-1,  -tstart-diff, -tstart, 10]), np.array([0, 0, Fuel_Ne, Fuel_Ne]))
    #plt.plot(t_fluxD, F_Ne)
    #plt.show()

    ss.eqsys.n_i.setFueling('Ne', F_Ne, times=t_fluxD)
    ss.eqsys.n_i.setFueling('D', fluxD*1.5e-1, times=t_fluxD) # ?


    # Radial grid
    ss.radialgrid.setB0(Btor)
    ss.radialgrid.setMinorRadius(a, t_Vp)
    ss.radialgrid.setMajorRadius(R0)
    ss.radialgrid.setWallRadius(0.4)
    ss.radialgrid.setVesselVolume(V_vessel)
    ss.radialgrid.setIref(2*np.pi*0.27*2e-3/scipy.constants.mu_0) # Maximal plasma current 1.2 MA, KSTAR 2 MA so same I_ref? Half? B is half for TCV ??
    ss.radialgrid.setBv(2.0e-3) # ??
    ss.radialgrid.setConnectionLengthFactor(connectionLengthFactor) # ??

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

    ss.other.include('fluid', 'stream', 'scalar')#, 'fluid/Tcold_ECH')

    #plt.plot(np.zeros(1), np.zeros(1))
    #plt.show()
    #ss.solver.setVerbose(True)
    #ss.solver.setDebug(savejacobian=True, savenumericaljacobian=True, timestep=2, iteration=5)

    return ss

def fun_to_min(c, stream, data):
    return c * stream - data

def drawplot1(axs, so, toffset=0, showlabel=False, save=False, first=True):
    """
    Draw a plot with a output from the given STREAMOutput object.
    """


    t = so.grid.t[:] + toffset
    Ip = so.eqsys.I_p[:, 0]
    Ire = so.eqsys.j_re.current()[:]
    ne = so.eqsys.n_cold[:, 0]
    Te = so.eqsys.T_cold[:, 0]

    hf = h5py.File('TCV65108.h5', 'r')
    t_Ip_d = np.array(hf.get('Plasma current').get('x')) + 0.01
    Ip_d = np.array(hf.get('Plasma current').get('z'))
    t_H = np.array(hf.get('H-alpha').get('x'))+ 0.01
    Halpha = np.array(hf.get('H-alpha').get('z'))
    Halpha_interp = interp1d(t_H, Halpha, 'linear')
    Halpha_i = Halpha_interp(t)
    t_C = np.array(hf.get('C-III').get('x'))+ 0.01
    CIII_d = np.array(hf.get('C-III').get('z'))
    CIII_interp = interp1d(t_C, CIII_d, 'linear')
    CIII_ds = CIII_interp(t)

    t_Energy = np.array(hf.get('Total energy (DML)').get('x')) + 0.01
    energy = np.array(hf.get('Total energy (DML)').get('z'))
    t_rad = np.array(hf.get('Radiated power').get('x')) + 0.01
    P_rad = np.array(hf.get('Radiated power').get('z')) * 1e3
    hf.close()
    hf = h5py.File('firdens.h5', 'r')
    t_FIR = np.array(hf.get('FIR').get('x')) + 0.01
    FIR = np.array(hf.get('FIR').get('z'))
    hf.close()
    hf = h5py.File('65108.h5', 'r')
    t_kin = np.array(hf.get('Kinetic energy').get('x')) + 0.01
    kin = np.array(hf.get('Kinetic energy').get('z'))
    t_T = np.array(hf.get('Average (Thomson)').get('x')) + 0.01
    temp = np.array(hf.get('Average (Thomson)').get('z'))
    hf.close()


    nD0 = so.eqsys.n_i['D'][0][:]
    nD1 = so.eqsys.n_i['D'][1][:]
    totD = nD0[1:].flatten() * so.other.stream.V_n_tot['D'][:].flatten() + nD1[1:].flatten() * so.other.stream.V_p[:,0].flatten()
    gammaD = (1 - nD0[1:].flatten() * so.other.stream.V_n_tot['D'][:].flatten() / totD)

    PEC_D = np.loadtxt('adas_H6561-9A.dat', unpack=True)
    PEC_D_ne = np.loadtxt('adas_Hne.dat', unpack=True) * 1e6
    PEC_D_Te = np.loadtxt('adas_HTe.dat', unpack=True)
    PEC_D_interp = interp2d(PEC_D_ne, PEC_D_Te, PEC_D)

    PEC_D_data = np.zeros(len(ne))
    for i in range(len(ne)):
        PEC_D_data[i] = PEC_D_interp(ne[i] , Te[i])
    Dalpha = ne * nD0.flatten() * PEC_D_data
    Dalpha *= np.max(Halpha)  / np.max(Dalpha)

    resD = least_squares(fun_to_min, 1, args=(Dalpha, Halpha_i))
    c_D = resD.x

    Dalpha *= c_D[0]

    PEC_C = np.loadtxt('adas_C4650-1.dat', unpack=True)
    PEC_C_ne = np.loadtxt('adas_Cne.dat', unpack=True) * 1e6
    PEC_C_Te = np.loadtxt('adas_CTe.dat', unpack=True)
    PEC_C_interp = interp2d(PEC_C_ne, PEC_C_Te, PEC_C)

    PEC_C_data = np.zeros(len(ne))
    for i in range(len(ne)):
        PEC_C_data[i] = PEC_C_interp(ne[i], Te[i])
    nC2 = so.eqsys.n_i['C'][2][:]# + so.eqsys.n_i['C'][4][:]
    CIII = ne * nC2.flatten() * PEC_C_data
    CIII *= np.max(CIII_ds) / np.max(CIII)

    resC = least_squares(fun_to_min, 1, args=(CIII, CIII_ds))
    c_C = resC.x

    CIII *= c_C[0]

    W_cold = so.eqsys.W_cold.integral()[:]
    P_cold = so.other.fluid.Tcold_radiation.integral()[:]

    Z_eff_num = nD1 * 1
    Z_eff_denom = nD1 * 1
    for i in range(1, 7):
        Z_eff_num += so.eqsys.n_i['C'][i][:] * i ** 2
        Z_eff_denom += so.eqsys.n_i['C'][i][:] * i
    Z_eff = Z_eff_num / Z_eff_denom

    eta_x = so.other.stream.eta_x[:]
    f_absorbedECH = 1 - np.exp(-eta_x / np.cos(10*np.pi/180))

    plotInternal(axs[0, 0], t_Ip_d, Ip_d / 1e6, ylabel=r'$I_{\rm p}$ (MA)', color='tab:blue', showlabel=showlabel, label='65108')
    plotInternal(axs[0, 0], t, Ip / 1e6, ylabel=r'$I_{\rm p}$ (MA)', color='tab:red', linestyle='dashed', showlabel=showlabel, label='STREAM')
    plotInternal(axs[0, 0], t, Ire / 1e6, ylabel=r'$I_{\rm p}$ (MA)', color='k', linestyle='dotted', showlabel=showlabel, label=r'$I_{\rm RE}$ STREAM')
    plotInternal(axs[0, 1], t_FIR, FIR / 1e19, ylabel=r'$n_{\rm e}$ ($1\cdot 10^{19}$m$^{-3}$)', color='tab:blue', showlabel=False, label='65108')
    plotInternal(axs[0, 1], t, ne / 1e19, ylabel=r'$n_{\rm e}$ ($1\cdot 10^{19}$m$^{-3}$)', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')
    plotInternal(axs[0, 2], t_T, temp, ylabel=r'$T_{\rm e}$ (eV)', color='tab:blue', showlabel=False, label='65108')
    plotInternal(axs[0, 2], t, Te, ylabel=r'$T_{\rm e}$ (eV)', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')
    #plotInternal(axs[1,1], t[1:], gammaD, ylabel=r'$\gamma_{\rm D}$ (\%)', color='tab:red', linestyle='dashed', showlabel=showlabel, label=r'$\gamma_{\rm D}$')
    plotInternal(axs[1, 0], t, Halpha_i, ylabel=r'${\rm D}_\alpha$ (a.u.)', color='tab:blue', showlabel=False, label='65108')
    plotInternal(axs[1, 0], t, Dalpha, ylabel=r'${\rm D}_\alpha$ (a.u.)', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')
    plotInternal(axs[1, 1], t, CIII_ds, ylabel=r'${\rm C}_{\rm III}$', color='tab:blue', showlabel=False, label='65108')
    plotInternal(axs[1, 1], t, CIII, ylabel=r'${\rm C}_{\rm III}$ (a.u.)', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')
    plotInternal(axs[1, 2], t_rad, P_rad/1e3, ylabel=r'$P_{\rm rad}$ [kW]', color='tab:blue', showlabel=False, label='65108')
    plotInternal(axs[1, 2], t[1:], P_cold/1e3, ylabel=r'$P_{\rm rad}$ [kW]', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')
    plotInternal(axs[2, 0], t_kin, kin, ylabel=r'$W_{\rm cold}$', color='tab:blue', showlabel=False, label='65108')
    plotInternal(axs[2, 0], t, W_cold, ylabel=r'$W_{\rm cold}$ [J]', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')
    plotInternal(axs[2 ,1], t, Z_eff, ylabel=r'$Z_{\rm eff}$', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')
    plotInternal(axs[2, 2], t[1:], f_absorbedECH*100, ylabel=r'$f_{\rm ECH, abs}$ (\%)', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')


    for i in range(axs.shape[0]):
        for j in range(axs.shape[1]):
            axs[i, j].set_xlim([0, 0.4])
            axs[i, j].grid(True, color='gainsboro', linewidth=0.5)
    #'''
    axs[0, 0].set_ylim([0, 0.4])
    axs[0, 1].set_ylim([0, 3.5])
    axs[0, 2].set_ylim([0, 700])
    axs[1, 0].set_ylim([0, 3.0])
    axs[1, 1].set_ylim([0, 0.05])
    axs[1, 2].set_ylim([0, 150])
    axs[2, 0].set_ylim([0, 6000])
    axs[2, 1].set_ylim([1, 1.5])
    axs[2, 2].set_ylim([0, 120])


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
    fig1, axs1 = plt.subplots(3, 3, figsize=(14, 12))

    drawplot1(axs1, so1, toffset=toffset)
    drawplot1(axs1, so2, toffset=toffset + so1.grid.t[-1], showlabel=True, first=False)
    axs1[0,0].legend(frameon=False)

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
        ss2.timestep.setTmax(0.4 - ss1.timestep.tmax)
        ss2.timestep.setNumberOfSaveSteps(0)
        ss2.timestep.setNt(100000)
        ss2.save(f'settings2{ext}.h5')
        so2 = runiface(ss2, f'output2{ext}.h5', quiet=False)
    else:
        so1 = STREAMOutput(f'output1{ext}.h5', loadsettings=False)
        so2 = STREAMOutput(f'output2{ext}.h5', loadsettings=False)

    if settings.plot:
        makeplots(so1, so2, toffset=0)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


