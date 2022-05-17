#!/usr/bin/env python3
#
# Script to reproduce the JET simulation in section 3 of (Kim et al., NF 2020).
# In order for STREAM results to match those of DYON/BKD0/SCENPLINT it is
# necessary to change the connection length parameter used in STREAM by
# editing the code. Locate the variable 'connectionLengthFactor' in
# 'include/STREAM/Equations/ConfinementTime.hpp' and change its value to '1'.
# Then, re-compile STREAM and run this script.

import argparse
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
from scipy.interpolate import interp2d
from scipy.constants import c, e, m_e, mu_0
import sys

sys.path.append('../../py')

#from run import makeplots
import PlasmaParameters as Formulas
import DREAM.Settings.Equations.ColdElectronTemperature as Tcold
import DREAM.Settings.Solver as Solver
import DREAM.Settings.Equations.ElectricField as ElectricField_D
import DREAM.Settings.Atomics as Atomics
from STREAM import STREAMOutput, STREAMSettings, runiface
import STREAM.Settings.Equations.ElectricField as ElectricField
import STREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways


def generate(prefill=1e-5, gamma=2e-2, fractionO = 0.001, fractionC = 0, Ip=1e4, tmax=1e-4, nt=2000, selfconsistent = False):
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
    n0 = 3e18*0.2 # Initial total deuterium density
    print(str(n0))
    nD = n0 * np.array([[1-gamma], [gamma]])

    Btor = 2.4      # Toroidal magnetic field [T]
    R0 = 2.9       # Plasma major radius [m]
    r_wall = 1      # Distance between plasma centre and passive structure [m]
    #r_wall = 1.3
    V_vessel = 175  # Vacuum vessel volume

    a = 0.8 # 1.25  # Plasma minor radius [m]

    tVloop = np.array([0.0, 0.005, 0.1, 0.11, 1.0, 2.0, 4.0])
    Vloop  = np.array([1.5, 1.3, 1.0, 0.5, 0.2, 0.15, 0.1])*11
    #plt.plot(tVloop, Vloop/10)
    #plt.xlim([-1,4])
    #plt.ylim([0,2])
    #plt.show()

    c1 = 0.98
    c2 = 0
    c3 = 1

    Te0 = 10     # electron temperature [eV]
    Ti0 = 1 # ion temperature [eV]

    # Initial electric field
    sigma = Formulas.evaluateSpitzerConductivity(n=nD[1], T=Te0, Z=1)
    j0 = Ip / (a ** 2 * np.pi)  # 943?
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
        Lp = float(mu_0 * R0 * (np.log(8 * R0 / a) + 0.25 - 2))
        print(str(Lp))
        ss.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
        ss.eqsys.E_field.setInitialProfile(efield=E0)
        ss.eqsys.E_field.setInductances(Lp=Lp, Lwall=9.1e-6, M=2.49e-6, Rwall=7.5e-4)
        ss.eqsys.E_field.setCircuitVloop(Vloop, tVloop)

    # Electron temperature
    ss.eqsys.T_cold.setType(Tcold.TYPE_SELFCONSISTENT)
    ss.eqsys.T_cold.setInitialProfile(Te0)
    ss.eqsys.T_cold.setRecombinationRadiation(recombination=Tcold.RECOMBINATION_RADIATION_INCLUDED)

    # Ions
    ss.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=nD, r=np.array([0]), T=Ti0)
    f = 5e18
    ss.eqsys.n_i.setFueling('D', [0, 0, f, f, 0, 0], [0, 0.4999, 0.5, 1.5, 1.5001, 4])

    # Disable runaway
    ss.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
    ss.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)

    # Radial grid
    ss.radialgrid.setB0(Btor)
    ss.radialgrid.setMinorRadius(a)
    ss.radialgrid.setMajorRadius(R0)
    ss.radialgrid.setWallRadius(r_wall)
    ss.radialgrid.setVesselVolume(V_vessel)
    ss.radialgrid.setBv(1e-3)

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


def drawplot1(axs, so, toffset=0, showlabel=False, save=True, first=True):
    """
    Draw a plot with a output from the given STREAMOutput object.
    """
    t = so.grid.t[:] + toffset
    Ip = so.eqsys.I_p[:, 0]
    Ire = so.eqsys.j_re.current()[:]
    ne = so.eqsys.n_cold[:, 0]
    Te = so.eqsys.T_cold[:, 0]

    PEC_D = np.loadtxt('adas_H6561-9A.dat', unpack=True)
    PEC_D_ne = np.loadtxt('adas_Hne.dat', unpack=True) * 1e6
    PEC_D_Te = np.loadtxt('adas_HTe.dat', unpack=True)
    PEC_D_interp = interp2d(PEC_D_ne, PEC_D_Te, PEC_D)

    nD0 = so.eqsys.n_i['D'][0][:]
    PEC_D_data = np.zeros(len(ne))
    for i in range(len(ne)):
        PEC_D_data[i] = PEC_D_interp(ne[i], Te[i])
    Dalpha = ne * nD0.flatten() * PEC_D_data
    Dalpha *= 0.8 / np.max(Dalpha)

    Ac = 1.25 ** 2 * np.pi
    streaming = Ip/(Ac*e*ne) * np.sqrt(m_e / (e * Te))

    plotInternal(axs[0], t, Ip / 1e6, ylabel=r'$I_{\rm p}$ (MA)', color='tab:blue', showlabel=False, label='STREAM')
    plotInternal(axs[0], t, Ire / 1e6, ylabel=r'$I_{\rm p}$ (MA)', color='tab:red', showlabel=False, label='STREAM')
    plotInternal(axs[1], t, streaming , ylabel=r'$\xi$ ', color='tab:blue', showlabel=False, label='STREAM')
    plotInternal(axs[2], t, Dalpha, ylabel=r'$D_\alpha$ ', color='tab:blue', showlabel=False, label='STREAM')
    plotInternal(axs[3], t, ne/1e18, ylabel=r'$n_{\rm e}$ ($10^{18}$ m$^{-3}$)', color='tab:blue', showlabel=False, label='STREAM')
    plotInternal(axs[4], t, Te/1e3, ylabel=r'$T_{\rm e}$ (keV)', color='tab:blue', showlabel=False, label='STREAM')


    for i in range(axs.shape[0]):
        axs[i].set_xlim([-1, 4])
        axs[i].grid(True)

    axs[0].set_ylim([0, 2])
    axs[1].set_ylim([0, 0.5])
    axs[2].set_ylim([0, 1])
    axs[3].set_ylim([0, 20])
    axs[4].set_ylim([0, 20])


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
    fig1, axs1 = plt.subplots(5, 1, figsize=(7, 10))

    drawplot1(axs1, so1)
    drawplot1(axs1, so2, toffset=so1.grid.t[-1], showlabel=True, first=False)


    fig1.tight_layout()
    #fig1.savefig('../../runaway_svn/mathias/startup/figs/deVries2020_f0.3.pdf')
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
        prefill = 2 * 2.7e-3 / 133.32    # Pa -> Torr
        ss1 = generate(prefill=prefill)
        ss1.save(f'settings1{ext}.h5')
        so1 = runiface(ss1, f'output1{ext}.h5', quiet=False)

        ss2 = STREAMSettings(ss1)
        ss2.fromOutput(f'output1{ext}.h5')
        ss2.timestep.setTmax(4 - ss1.timestep.tmax)
        ss2.timestep.setNt(50000)
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


