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

def generate(nD0=5e17, gamma=2e-3, P_inj=665e3, theta=10, phi=90, N=2, f_x=1.0, f_o=0.0, tmax=1e-5, nt=10000, tstart=0.0, texpstart=0.0, inputdatafilename='', C_Vloop=1.0, C_flux=1.5e-1):
    """
    Generate a STREAMSettings object for a simulation with the specified
    parameters.

    :param nD0:                Total initial deuterium density [1/m^3]
    :param gamma:              Degree of ionization
    :param P_inj:              Injected power [W]
    :param f_x:                X mode fraction
    :param f_o:                O mode fraction
    :param theta:              ECH poloidal launch angle [deg]
    :param phi:                ECH toroidal launch angle [deg]
    :param N:                  Fundamental harmonic
    :param tmax:               Simulation time of first simulation [s]
    :param nt:                 Number of time steps
    :param tstart:             Simulation start time in TCV time data [s]
    :param texpstart:          Experiment start time in TCV time data [s]
    :param inputdatafilename:  Name of TCV input data file
    :param C_Vloop:            Loop voltage factor
    :param C_flux:             Deuterium flux factor
    """
    if inputdatafilename == '':
        raise Exception('Filename for TCV data is missing')
    elif len(inputdatafilename) <= 3 or inputdatafilename[-3:] != '.h5':
        inputdatafilename += '.h5'

    hf = h5py.File(inputdatafilename, 'r')

    # Initial total deuterium density
    nD = nD0 * np.array([[1 - gamma], [gamma]])  # Deuterium density [1/m^3]
    nC = 1e3                                     # Carbon density [1/m^3]

    # Initial temperatures
    Te0 = 1                                      # Electron temperature [eV]
    Ti0 = 0.026                                  # Ion temperature [eV]

    # ECH parameters
    theta *= np.pi / 180                         # ECH poloidal launch angle [rad]
    phi   *= np.pi / 180                         # ECH toroidal launch angle [rad]

    # Tokamak geometric parameters
    R0 = 0.89                                    # Plasma major radius [m]
    Btor = 1.45                                  # Toroidal magnetic field [T]
    V_vessel = 4.632                             # Vacuum vessel volume

    t_Vp_osc = np.array(hf.get('Plasma volume').get('x'))                                # TCV time data for plasma volume [s]
    V_p_osc = np.array(hf.get('Plasma volume').get('z'))                                 # TCV plasma volume data [m^3]
    V_p_s = savgol_filter(V_p_osc, 65, 3)                                                # Smoothed plasma volume [m^3]
    V_initfun = interp1d(np.append(np.array([texpstart - tstart]), t_Vp_osc - tstart),
                         np.append(np.array([V_vessel]), V_p_s), 'cubic')                #
    t_a = np.linspace(0, t_Vp_osc[-1] - tstart, 200)                                    # Shifted minor radius time vector [s]
    V_p = V_initfun(t_a)                                                                # Smoothed plasma volume data starting at vessel volume [m^3]
    a = np.sqrt(V_p / (2 * np.pi ** 2 * R0))                                             # Plasma minor radius [m]

    # Loop voltage (smoothed TCV data)
    t_loop     = np.array(hf.get('Loop voltage').get('x'))                                 # TCV time data for loop voltage [s]
    V_loop_osc =- np.array(hf.get('Loop voltage').get('z'))                                # TCV loop voltage data [V]
    it0_loop   = int((tstart - t_loop[0]) / (t_loop[-1] - t_loop[0]) * t_loop.shape[0])-1  # Index of tstart in TCV time data for loop voltage
    t_loop     = t_loop[it0_loop:]-tstart                                                  # Shifted loop voltage time vector [s]
    V_loop     = savgol_filter(V_loop_osc[it0_loop:], 157, 3)*C_Vloop                      # Smoothed loop voltage [V]
    plt.plot(t_loop, V_loop)
    plt.show()

    # Initial electric field (derived from TCV data)
    t_Ip = np.array(hf.get('Plasma current').get('x'))                           # TCV time data for plasma current [s]
    Ip = np.array(hf.get('Plasma current').get('z'))                             # TCV plasma current data [A]
    it0_Ip = int((tstart - t_Ip[0]) / (t_Ip[-1] - t_Ip[0]) * t_Ip.shape[0]) - 1  # Index of tstart in TCV time data for plasma current
    Ip = Ip[it0_Ip]                                                              # Initial plasma current [A]
    j0 = Ip / (a[0] ** 2 * np.pi)                                                # Initial plasma current density [A/m^2]
    E0 = j0 / Formulas.evaluateSpitzerConductivity(n=nD[1], T=Te0, Z=1)          # Initial electric field [V/m^2]

    # Deuterium fueling data
    t_fluxD = np.array(hf.get('Particle flux (D2)').get('x'))                                # TCV time data for deuterium flux [s]
    fluxD = np.array(hf.get('Particle flux (D2)').get('z'))                                  # TCV deuterium flux data [1/m^3 s]
    it0_flux = int((tstart - t_fluxD[0]) / (t_fluxD[-1] - t_fluxD[0]) * t_fluxD.shape[0])-1  # Index of tstart in TCV deuterium flux data for plasma current
    t_fuel = t_fluxD[it0_flux:]-tstart                                                       # Shifted deuterium fueling time vector [s]
    fuel = fluxD[it0_flux:]*C_flux                                                           # Deuterium fueling vector [1/m^3 s]

    # Other parameters
    l_i  = 0.0  # Internal inductance
    C_Lf = 3.0  # Connection lengh factor

    hf.close()

    ss = STREAMSettings()

    ss.atomic.adas_interpolation = Atomics.ADAS_INTERP_BILINEAR

    # Electric field
    ss.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
    ss.eqsys.E_field.setInitialProfile(E0)
    Lp = float(scipy.constants.mu_0 * R0 * (np.log(8 * R0 / 0.25) + l_i/2 - 2))
    ss.eqsys.E_field.setInductances(Lp=Lp, Lwall=9.1e-6, M=2.49e-6, Rwall=1e6)
    ss.eqsys.E_field.setCircuitVloop(V_loop, t_loop)

    # Electron temperature
    ss.eqsys.T_cold.setType(Tcold.TYPE_SELFCONSISTENT)
    ss.eqsys.T_cold.setInitialProfile(Te0)
    ss.eqsys.T_cold.setRecombinationRadiation(recombination=Tcold.RECOMBINATION_RADIATION_INCLUDED)

    # Ions
    ss.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=nD, r=np.array([0]), T=Ti0)
    ss.eqsys.n_i.addIon(name='C', Z=6, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=nC, r=np.array([0]), T=Ti0)

    # Disable runaway
    ss.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
    ss.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)

    # Recycling coefficients
    iD = ss.eqsys.n_i.getIndex('D')
    ss.eqsys.n_i.ions[iD].setRecyclingCoefficient('D', 1)

    ss.eqsys.n_i.ions[iD].setRecyclingCoefficient('C', [0.015, 0.05, 0.015, 0.010], [0, 0.00001, 0.05, 0.4])  # ?

    # Fueling
    ss.eqsys.n_i.setFueling('D', fuel, times=t_fuel) # ?

    # Radial grid
    ss.radialgrid.setB0(Btor)
    ss.radialgrid.setMinorRadius(a, t_a)
    ss.radialgrid.setMajorRadius(R0)
    ss.radialgrid.setWallRadius(0.4)
    ss.radialgrid.setVesselVolume(V_vessel)
    ss.radialgrid.setIref(2*np.pi*0.27*2e-3/scipy.constants.mu_0)
    print(2*np.pi*0.27*2e-3/scipy.constants.mu_0)
    ss.radialgrid.setBv(2.0e-3)
    ss.radialgrid.setConnectionLengthFactor(C_Lf)

    # ECH
    ss.radialgrid.setECHParameters(P_inj, f_o, f_x, theta, phi, N)

    # Disable kinetic grids
    ss.hottailgrid.setEnabled(False)
    ss.runawaygrid.setEnabled(False)

    # Numerical settings
    ss.solver.setType(Solver.NONLINEAR)
    ss.solver.preconditioner.setEnabled(False)
    ss.timestep.setTmax(tmax)
    ss.timestep.setNt(nt)

    ss.other.include('fluid', 'stream', 'scalar', 'fluid/Tcold_ECH')

    return ss

def fun_to_min(c, stream, data):
    return c * stream - data

def drawplotSimulation(axs, so, toffset=0, tstart=0, tend=1.0, showlabel=False, filename=''):
    """
    Draw a plot with an output from the given STREAMOutput object.
    """
    # STREAM output parameters
    t = so.grid.t[:] + toffset
    Ip = so.eqsys.I_p[:, 0]
    Ire = so.eqsys.j_re.current()[:]
    ne = so.eqsys.n_cold[:, 0]
    Te = so.eqsys.T_cold[:, 0]
    nD0 = so.eqsys.n_i['D'][0][:]
    nD1 = so.eqsys.n_i['D'][1][:]
    nC2 = so.eqsys.n_i['C'][2][:]
    W_cold = so.eqsys.W_cold.integral()[:]
    P_cold = so.other.fluid.Tcold_radiation.integral()[:]

    hf = h5py.File(filename, 'r')

    # Deuterium-alpha emmission evaluation
    t_H = np.array(hf.get('H-alpha').get('x'))-tstart
    Halpha = np.array(hf.get('H-alpha').get('z'))
    Halpha_interp = interp1d(t_H, Halpha, 'linear')
    Halpha_i = Halpha_interp(t)
    PEC_D = np.loadtxt('adas_H6561-9A.dat', unpack=True)
    PEC_D_ne = np.loadtxt('adas_Hne.dat', unpack=True) * 1e6
    PEC_D_Te = np.loadtxt('adas_HTe.dat', unpack=True)
    PEC_D_interp = interp2d(PEC_D_ne, PEC_D_Te, PEC_D)
    PEC_D_data = np.zeros(len(ne))
    for i in range(len(ne)):
        PEC_D_data[i] = PEC_D_interp(ne[i], Te[i])
    Dalpha = ne * nD0.flatten() * PEC_D_data
    Dalpha *= np.max(Halpha) / np.max(Dalpha)
    resD = least_squares(fun_to_min, 1, args=(Dalpha, Halpha_i))
    c_D = resD.x
    Dalpha *= c_D[0]

    # Carbon-III emmission evaluation
    t_C = np.array(hf.get('C-III').get('x'))-tstart
    CIII_d = np.array(hf.get('C-III').get('z'))
    CIII_interp = interp1d(t_C, CIII_d, 'linear')
    CIII_ds = CIII_interp(t)
    hf.close()
    PEC_C = np.loadtxt('adas_C4650-1.dat', unpack=True)
    PEC_C_ne = np.loadtxt('adas_Cne.dat', unpack=True) * 1e6
    PEC_C_Te = np.loadtxt('adas_CTe.dat', unpack=True)
    PEC_C_interp = interp2d(PEC_C_ne, PEC_C_Te, PEC_C)
    PEC_C_data = np.zeros(len(ne))
    for i in range(len(ne)):
        PEC_C_data[i] = PEC_C_interp(ne[i], Te[i])
    CIII = ne * nC2.flatten() * PEC_C_data
    CIII *= np.max(CIII_ds) / np.max(CIII)
    resC = least_squares(fun_to_min, 1, args=(CIII, CIII_ds))
    c_C = resC.x
    CIII *= c_C[0]

    plotInternal(axs[0], t + tstart, Ip / 1e6, ylabel=r'$I_{\rm p}$ (MA)', color='tab:red', linestyle='dashed', showlabel=showlabel, label='STREAM')
    plotInternal(axs[0], t + tstart, Ire / 1e6, ylabel=r'$I_{\rm p}$ (MA)', color='k', linestyle='dotted', showlabel=showlabel, label=r'$I_{\rm RE}$ STREAM')
    plotInternal(axs[1], t + tstart, Te, ylabel=r'$T_{\rm e}$ (eV)', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')
    plotInternal(axs[2], t + tstart, ne / 1e19, ylabel=r'$n_{\rm e}$ ($1\cdot 10^{19}$m$^{-3}$)', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')
    plotInternal(axs[3], t[1:], P_cold / 1e3, ylabel=r'$P_{\rm rad}$ [kW]', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')
    plotInternal(axs[4], t + tstart, Dalpha, ylabel=r'${\rm D}_\alpha$ (a.u.)', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')
    plotInternal(axs[5], t + tstart, CIII, ylabel=r'${\rm C}_{\rm III}$ (a.u.)', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')
    plotInternal(axs[6], t + tstart, W_cold, ylabel=r'$W_{\rm cold}$ [J]', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')

    for i in range(len(axs)):
        axs[i].set_xlim([tstart, tend])
        axs[i].grid(True, color='gainsboro', linewidth=0.5)


def drawplotTCV(axs, tstart=0, tend=1.0, showlabel=False, filename=''):
    """
    Draw a plot with an output from the given STREAMOutput object.
    """
    hf = h5py.File(filename, 'r')
    t_Ip_d = np.array(hf.get('Plasma current').get('x'))
    Ip_d = np.array(hf.get('Plasma current').get('z'))
    t_H = np.array(hf.get('H-alpha').get('x'))
    Halpha = np.array(hf.get('H-alpha').get('z'))
    t_C = np.array(hf.get('C-III').get('x'))
    CIII_d = np.array(hf.get('C-III').get('z'))
    t_Energy = np.array(hf.get('Total energy (DML)').get('x'))
    energy = np.array(hf.get('Total energy (DML)').get('z'))
    t_rad = np.array(hf.get('Radiated power').get('x'))
    P_rad = np.array(hf.get('Radiated power').get('z')) * 1e3
    hf.close()
    hf = h5py.File('firdens.h5', 'r')                         # TODO: Fix filename
    t_FIR = np.array(hf.get('FIR').get('x'))
    FIR = np.array(hf.get('FIR').get('z'))
    hf.close()
    hf = h5py.File('65108.h5', 'r')                           # TODO: Fix filename
    t_kin = np.array(hf.get('Kinetic energy').get('x'))
    kin = np.array(hf.get('Kinetic energy').get('z'))
    t_T = np.array(hf.get('Average (Thomson)').get('x'))
    temp = np.array(hf.get('Average (Thomson)').get('z'))
    hf.close()

    plotInternal(axs[0], t_Ip_d, Ip_d / 1e6, ylabel=r'$I_{\rm p}$ (MA)', color='tab:blue', showlabel=showlabel, label='65108')
    plotInternal(axs[1], t_T, temp, ylabel=r'$T_{\rm e}$ (eV)', color='tab:blue', showlabel=False, label='65108')
    plotInternal(axs[2], t_FIR, FIR / 1e19, ylabel=r'$n_{\rm e}$ ($1\cdot 10^{19}$m$^{-3}$)', color='tab:blue', showlabel=False, label='65108')
    plotInternal(axs[3], t_rad, P_rad / 1e3, ylabel=r'$P_{\rm rad}$ [kW]', color='tab:blue', showlabel=False, label='65108')
    plotInternal(axs[4], t_H, Halpha, ylabel=r'${\rm D}_\alpha$ (a.u.)', color='tab:blue', showlabel=False, label='65108')
    plotInternal(axs[5], t_C, CIII_d, ylabel=r'${\rm C}_{\rm III}$', color='tab:blue', showlabel=False, label='65108')
    plotInternal(axs[6], t_kin, kin, ylabel=r'$W_{\rm cold}$', color='tab:blue', showlabel=False, label='65108')

    for i in range(len(axs)):
        axs[i].set_xlim([tstart, tend])
        axs[i].grid(True, color='gainsboro', linewidth=0.5)


def drawplotOther(axs, so, toffset=0, tstart=0, tend=1.0, showlabel=False, theta=10):
    """
    Draw a plot with an output from the given STREAMOutput object.
    """
    # STREAM output parameters
    t = so.grid.t[:] + toffset
    nD0 = so.eqsys.n_i['D'][0][:]
    nD1 = so.eqsys.n_i['D'][1][:]

    # Absorbed ECH power
    eta_x = so.other.stream.eta_x[:]
    f_absorbedECH = 1 - np.exp(-eta_x / np.cos(theta*np.pi/180))

    # Effective charge
    Z_eff_num = nD1 * 1
    Z_eff_denom = nD1 * 1
    for i in range(1, 7):
        Z_eff_num   += so.eqsys.n_i['C'][i][:] * i ** 2
        Z_eff_denom += so.eqsys.n_i['C'][i][:] * i
    Z_eff = Z_eff_num / Z_eff_denom

    # Degree of ionization
    totD = nD0[1:].flatten() * so.other.stream.V_n_tot['D'][:].flatten() + nD1[1:].flatten() * so.other.stream.V_p[:,0].flatten()
    gammaD = (1 - nD0[1:].flatten() * so.other.stream.V_n_tot['D'][:].flatten() / totD)

    plotInternal(axs[0], t[1:] + tstart, f_absorbedECH * 100, ylabel=r'$f_{\rm ECH, abs}$ (\%)', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')
    plotInternal(axs[1], t + tstart, Z_eff, ylabel=r'$Z_{\rm eff}$', color='tab:red', linestyle='dashed', showlabel=False, label='STREAM')
    plotInternal(axs[2], t[1:] + tstart, gammaD*100, ylabel=r'$\gamma_{\rm D}$ (\%)', color='tab:red', linestyle='dashed', showlabel=showlabel, label=r'$\gamma_{\rm D}$')


    for i in range(axs.shape[0]):
        axs[i].set_xlim([tstart, tend])
        axs[i].grid(True, color='gainsboro', linewidth=0.5)

    axs[0].set_ylim([0, 110])
    axs[1].set_ylim([0, 2])
    axs[2].set_ylim([0, 110])

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


def makeplots(so1, so2, filename='', tstart=0.0, tend=1.0, theta=10, savefig=False, pathfig=''):
    if pathfig != '' and pathfig[-1] != '/':
        pathfig += '/'
    if filename == '':
        raise Exception('Filename for TCV data is missing')
    elif len(filename) <= 3 or filename[-3:] != '.h5':
        if savefig:
            fig1filename = pathfig + filename + '_vs.pdf'
            fig2filename = pathfig + filename + '_other.pdf'
        filename += '.h5'
    elif savefig:
        fig1filename = pathfig + filename[:-3] + '_vs.pdf'
        fig2filename = pathfig + filename[:-3] + '_other.pdf'

    fig1 = plt.figure(figsize=(14, 12))
    axs1 = [fig1.add_subplot(3, 3, (1, 2)), fig1.add_subplot(3, 3, 3),
            fig1.add_subplot(3, 3, (4, 5)), fig1.add_subplot(3, 3, 6),
            fig1.add_subplot(3, 3, 7), fig1.add_subplot(3, 3, 8), fig1.add_subplot(3, 3, 9)]

    drawplotTCV(axs1, tstart=tstart, tend=tend, showlabel=True, filename=filename)
    drawplotSimulation(axs1, so1, toffset=0, tstart=tstart, tend=tend, filename=filename)
    drawplotSimulation(axs1, so2, toffset=so1.grid.t[-1], tstart=tstart, tend=tend, showlabel=True, filename=filename)
    axs1[0].legend(frameon=False)
    for i in range(len(axs1)):
        axs1[i].set_ylim(bottom=0.0)
    axs1[4].set_ylim(top=2.7)
    fig1.tight_layout()
    if savefig:
        fig1.savefig(fig1filename)

    fig2, axs2 = plt.subplots(1, 3, figsize=(14, 4))
    drawplotOther(axs2, so1, toffset=0, tstart=tstart, tend=tend, theta=theta)
    drawplotOther(axs2, so2, toffset=so1.grid.t[-1], tstart=tstart, tend=tend, theta=theta)
    fig2.tight_layout()
    if savefig:
        fig2.savefig(fig2filename)
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

    tstart = -0.01
    tend   = 0.4
    if not settings.skip:
        ss1 = generate(nD0 = 5e17, tmax = 1e-5, nt = 10000, tstart = -0.01, texpstart = -0.02, inputdatafilename = 'TCV65108_input.h5', C_Vloop = 1.0, C_flux = 1.5e-1)
        ss1.save(f'settings1{ext}.h5')
        so1 = runiface(ss1, f'output1{ext}.h5', quiet=False)

        ss2 = STREAMSettings(ss1)
        ss2.fromOutput(f'output1{ext}.h5')
        ss2.timestep.setTmax(tend - tstart - ss1.timestep.tmax)
        ss2.timestep.setNumberOfSaveSteps(0)
        ss2.timestep.setNt(30000)
        ss2.save(f'settings2{ext}.h5')
        so2 = runiface(ss2, f'output2{ext}.h5', quiet=False)
    else:
        so1 = STREAMOutput(f'output1{ext}.h5', loadsettings=False)
        so2 = STREAMOutput(f'output2{ext}.h5', loadsettings=False)

    if settings.plot:
        makeplots(so1, so2, filename='TCV65108.h5', tstart=tstart, tend=tend)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


