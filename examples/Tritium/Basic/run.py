#!/usr/bin/env python3
#
# This example uses parameters from Table 4.1 of [H. T. Kim, PhD Thesis (2013)]
# to simulate basic Deuterium burn-through in a JET-like plasma.
#

import argparse
import matplotlib.pyplot as plt
import numpy as np
import sys
import scipy.constants
sys.path.append('../../../py')

import PlasmaParameters as Formulas
import DREAM.Settings.Equations.ColdElectronTemperature as Tcold
import DREAM.Settings.Solver as Solver
import DREAM.Settings.Equations.ElectricField as ElectricField_D
import DREAM.Settings.Atomics as Atomics
from STREAM import STREAMOutput, STREAMSettings, runiface
import STREAM.Settings.Equations.ElectricField as ElectricField
import STREAM.Settings.Equations.IonSpecies as Ions


def generate(prefill=5e-5, gamma=2e-3, fractionT=0, Vloop=12, Vloop_t=0, j0=883.3, tmax=1e-4, nt=4000, EfieldDYON = True):
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
    nD = n0 * (1-fractionT) * np.array([[1-gamma], [gamma]])
    nT = n0 * fractionT * np.array([[1-gamma], [gamma]])

    Btor = 2.65     # Toroidal magnetic field [T]
    a = 1.6         # Plasma minor radius [m]
    R0 = 5.65       # Plasma major radius [m]
    l_MK2 = 1       # Distance between plasma centre and passive structure [m]
    V_vessel = 100  # Vacuum vessel volume

    Te0 = 1     # electron temperature [eV]
    Ti0 = 0.026 # ion temperature [eV]

    # Initial electric field
    E0 = j0 / Formulas.evaluateSpitzerConductivity(n=nD[1], T=Te0, Z=1)

    ss = STREAMSettings()

    ss.atomic.adas_interpolation = Atomics.ADAS_INTERP_BILINEAR

    # Electric field
    if EfieldDYON:
        Lp = float(scipy.constants.mu_0 * R0 * (np.log(8 * R0 / a) + 0.25 - 2))

        ss.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
        ss.eqsys.E_field.setInitialProfile(E0)
        ss.eqsys.E_field.setInductances(Lp=Lp, Lwall=9.1e-6, M=2.49e-6, Rwall=1e6)
        ss.eqsys.E_field.setCircuitVloop(Vloop, Vloop_t)
    else:
        R = 1e6
        L = 9.1e-6
        wall_time = L / R

        ss.eqsys.E_field.setType(ElectricField_D.TYPE_SELFCONSISTENT)
        ss.eqsys.E_field.setInitialProfile(efield=E0)
        ss.eqsys.E_field.setBoundaryCondition(ElectricField_D.BC_TYPE_TRANSFORMER, V_loop_wall_R0=Vloop/R0, times=Vloop_t, inverse_wall_time=1/wall_time, R0=R0)

    # Electron temperature
    ss.eqsys.T_cold.setType(Tcold.TYPE_SELFCONSISTENT)
    ss.eqsys.T_cold.setInitialProfile(Te0)
    ss.eqsys.T_cold.setRecombinationRadiation(recombination=Tcold.RECOMBINATION_RADIATION_INCLUDED)

    # Ions
    ss.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=nD, r=np.array([0]), T=Ti0)
    if fractionT > 0:
        ss.eqsys.n_i.addIon(name='T', Z=1, iontype=Ions.IONS_DYNAMIC, n=nT, r=np.array([0]), T=Ti0, tritium=True)

    # Disable runaway
    ss.eqsys.n_re.setAvalanche(False)
    ss.eqsys.n_re.setDreicer(False)

    # Recycling coefficients
    iD = ss.eqsys.n_i.getIndex('D')
    ss.eqsys.n_i.ions[iD].setRecyclingCoefficient('D', 1)

    iT = ss.eqsys.n_i.getIndex('T')
    ss.eqsys.n_i.ions[iT].setRecyclingCoefficient('T', 1)

    # Radial grid
    ss.radialgrid.setB0(Btor)
    ss.radialgrid.setMinorRadius(a)
    ss.radialgrid.setMajorRadius(R0)
    ss.radialgrid.setWallRadius(l_MK2)
    ss.radialgrid.setVesselVolume(V_vessel)
    
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


def drawplot1(axs, so, color='r', toffset=0):
    """
    Draw a plot with a output from the given STREAMOutput object.
    """
    t = so.grid.t[:] + toffset
    plotInternal(axs[0], t[1:], np.diff(so.eqsys.W_cold[:,0]) / np.diff(so.grid.t[:]), ylabel=r'Power consumption (W/m$^3$)', color=color, ylim=[0,50])
    plotInternal(axs[1], t, so.eqsys.I_p[:,0], ylabel=r'Plasma current $I_{\rm p}$ (A)', color=color, log=True)


def drawplot2(axs, so, color='r', toffset=0):
    t = so.grid.t[:] + toffset

    V_p = so.other.stream.V_p[:,0]
    V_n_tot = so.other.stream.V_n_tot['D'][:]
    #dWr_dt  = np.diff(so.other.fluid.Tcold_radiation[:,0]) / np.diff(t[1:])
    nD0 = so.eqsys.n_i['D'][0][1:,0]
    nD1 = so.eqsys.n_i['D'][1][1:,0]
    #gamma_i = nD1*V_p / (nD1*V_p + nD0*V_n_tot)
    gamma_i = nD1 / (nD1 + nD0)

    plotInternal(axs[0], t[1:], so.other.fluid.Tcold_radiation[:,0], 'Power loss (W)', color=color, xlbl=False)
    plotInternal(axs[1], t[1:], gamma_i*100, r'Degree of ionization (\%)', color=color, ylim=[0,105], xlbl=False)
    plotInternal(axs[2], t, so.eqsys.T_cold[:,0], 'Electron temperature (eV)', color=color, xlbl=False)
    plotInternal(axs[3], t, so.eqsys.n_cold[:,0]*1e-18, 'Electron density (10$^{18}$ m$^-3$)', color=color)


def drawplot3(axs, so, toffset=0, showlabel=True):
    t = so.grid.t[:] + toffset

    V_p = so.other.stream.V_p[:,0]

    Prad    = so.other.fluid.Tcold_radiation[:,0] * V_p
    Pequi   = so.other.fluid.Tcold_ion_coll[:,0] * V_p
    Ptransp = so.other.scalar.energyloss_T_cold[:,0] * V_p
    P_tot   = Prad + Pequi + Ptransp
    P_net   = P_tot + so.other.fluid.Tcold_ohmic[:,0] * V_p

    nD0 = so.eqsys.n_i['D'][0][1:,0]
    nD1 = so.eqsys.n_i['D'][1][1:,0]
    gamma_i = nD1 / (nD1 + nD0)

    plotInternal(axs[0], t, so.eqsys.I_p[:,0]/1e3, r'$I_{\rm p}$ (kA)', linestyle='-.', color='b', xlbl=False)
    plotInternal(axs[1], t[1:], gamma_i*100, r'Degree of ionization (\%)', color='k', ylim=[0,105], xlbl=False)

    ylbl = r'Power balance'
    plotInternal(axs[2], t[1:], P_tot, ylbl, label='Total electron power loss', showlabel=showlabel, color='r', linestyle='--')
    plotInternal(axs[2], t[1:], Prad, ylbl, label='Radiation + ionization', showlabel=showlabel, color='b', linestyle='-.')
    plotInternal(axs[2], t[1:], Pequi, ylbl, label='Equilibration', showlabel=showlabel, color='g', linestyle='--')
    plotInternal(axs[2], t[1:], Ptransp, ylbl, label='Electron transport', showlabel=showlabel, color='m', linestyle=':')
    plotInternal(axs[2], t[1:], P_net, ylbl, label='Net electron heating power', showlabel=showlabel, color='k', linestyle='-', ylim=[0, 1.5e5])

    axs[0].set_xlim([0, 0.03])
    axs[1].set_xlim([0, 0.03])
    axs[2].set_xlim([0, 0.03])
    axs[2].legend(frameon=False)


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


def makeplots(so11, so12, so21, so22):
    fig1, axs1 = plt.subplots(2, 1, figsize=(7,5), sharex=True)

    drawplot1(axs1, so11, color='b')
    drawplot1(axs1, so12, color='b', toffset=so11.grid.t[-1])
    drawplot1(axs1, so21, color='r')
    drawplot1(axs1, so22, color='r', toffset=so21.grid.t[-1])
    axs1[0].legend([r'D', '', r'DT'])

    fig2, axs2 = plt.subplots(4, 1, figsize=(7, 10))
    drawplot2(axs2, so11, color='b')
    drawplot2(axs2, so12, color='b', toffset=so11.grid.t[-1])
    drawplot2(axs2, so21, color='r')
    drawplot2(axs2, so22, color='r', toffset=so21.grid.t[-1])

    axs2[0].legend([r'D', '', r'DT'])

    fig3, axs3 = plt.subplots(3, 1, figsize=(7, 10))
    fig4, axs4 = plt.subplots(3, 1, figsize=(7, 10))

    drawplot3(axs3, so11)
    drawplot3(axs3, so12, toffset=so11.grid.t[-1], showlabel=False)
    fig3.suptitle('D')
    drawplot3(axs4, so21)
    drawplot3(axs4, so22, toffset=so21.grid.t[-1], showlabel=False)
    fig4.suptitle('DT')

    plt.tight_layout()
    plt.show()


def main(argv):
    FONTSIZE = 16
    plt.rcParams.update({'font.size': FONTSIZE})

    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--skip-all', help="Skip all simulations and only generate plots", dest="skip", action='store_const', const=[1,2])
    parser.add_argument('-e', '--extension', help="Append the specified string to the end of file names", dest="extension", action='store', default='')
    parser.add_argument('-n', '--no-plot', help="Only run simulations, don't generate plots", dest="plot", action='store_false', default=True)
    parser.add_argument('-s', '--skip', help="Skip the specified simulation(s)", dest="skip", nargs="*", type=int)

    settings = parser.parse_args()

    ext = '' if not settings.extension else '_' + settings.extension

    if settings.skip is None or (len(settings.skip) > 0 and 1 not in settings.skip):
        print('RUN 1')
        ss11 = generate(EfieldDYON=False)
        ss11.save(f'settings1WithoutT{ext}.h5')
        so11 = runiface(ss11, f'output1WithoutT{ext}.h5', quiet=False)

        ss12 = STREAMSettings(ss11)
        ss12.fromOutput(f'output1WithoutT{ext}.h5')
        ss12.timestep.setTmax(0.03 - ss11.timestep.tmax)
        ss12.timestep.setNt(3000)
        ss12.save(f'settings2WithoutT{ext}.h5')
        so12 = runiface(ss12, f'output2WithoutT{ext}.h5', quiet=False)
    else:
        so11 = STREAMOutput(f'output1WithoutT{ext}.h5')
        so12 = STREAMOutput(f'output2WithoutT{ext}.h5')

    if settings.skip is None or (len(settings.skip) > 0 and 2 not in settings.skip):
        print('RUN 2')
        ss21 = generate(fractionT=0.5, EfieldDYON=False)
        ss21.save(f'settings1WithT{ext}.h5')
        so21 = runiface(ss21, f'output1WithT{ext}.h5', quiet=False)

        ss22 = STREAMSettings(ss21)
        ss22.fromOutput(f'output1WithT{ext}.h5')
        ss22.timestep.setTmax(0.03 - ss21.timestep.tmax)
        ss22.timestep.setNt(3000)
        ss22.save(f'settings2WithT{ext}.h5')
        so22 = runiface(ss22, f'output2WithT{ext}.h5', quiet=False)
    else:
        so21 = STREAMOutput(f'output1WithT{ext}.h5')
        so22 = STREAMOutput(f'output2WithT{ext}.h5')

    if settings.plot:
        makeplots(so11, so12, so21, so22)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


