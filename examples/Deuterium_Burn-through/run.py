#!/usr/bin/env python3
#
# This example uses parameters from Table 4.1 of [H. T. Kim, PhD Thesis (2013)]
# to simulate basic Deuterium burn-through in a JET-like plasma.
#

import argparse
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../py')

import PlasmaParameters as Formulas
import DREAM.Settings.Equations.ColdElectronTemperature as Tcold
import DREAM.Settings.Solver as Solver

from STREAM import STREAMOutput, STREAMSettings, runiface
import STREAM.Settings.Equations.ElectricField as ElectricField
import STREAM.Settings.Equations.IonSpecies as Ions


def generate(prefill=5e-5, gamma=2e-3, Vloop=20, Vloop_t=0, j0=405.8, tmax=0.003, nt=1000):
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
    n0 = 3.22e22 * prefill  # Initial total deuterium density
    nD = n0 * np.array([[1-gamma], [gamma]])

    Btor = 2.3      # Toroidal magnetic field [T]
    a = 0.5         # Plasma minor radius [m]
    R0 = 3          # Plasma major radius [m]
    l_MK2 = 1       # Distance between plasma centre and passive structure [m]
    V_vessel = 100  # Vacuum vessel volume

    Te0 = 1     # electron temperature [eV]
    Ti0 = 0.03  # ion temperature [eV]

    # Initial electric field
    E0 = j0 / Formulas.evaluateSpitzerConductivity(n=nD[1], T=Te0, Z=1)

    ss = STREAMSettings()

    # Electric field
    ss.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
    ss.eqsys.E_field.setInitialProfile(E0)
    #ss.eqsys.E_field.setInductances(Lp=6.09e-6, Lwall=9.1e-6, M=2.49e-6, Rwall=1e6)
    ss.eqsys.E_field.setInductances(Lp=8e-6, Lwall=9.1e-6, M=2.49e-6, Rwall=1e6)
    ss.eqsys.E_field.setCircuitVloop(Vloop, Vloop_t)

    # Electron temperature
    ss.eqsys.T_cold.setType(Tcold.TYPE_SELFCONSISTENT)
    ss.eqsys.T_cold.setInitialProfile(Te0)
    ss.eqsys.T_cold.setRecombinationRadiation(recombination=Tcold.RECOMBINATION_RADIATION_INCLUDED)

    # Ions
    ss.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=nD, r=np.array([0]), T=Ti0)

    # Disable runaway
    ss.eqsys.n_re.setAvalanche(False)
    ss.eqsys.n_re.setDreicer(False)

    # Recycling coefficients (unused)
    ss.eqsys.n_i.setJET_CWrecycling()

    # Radial grid
    ss.radialgrid.setB0(Btor)
    ss.radialgrid.setMinorRadius(a)
    ss.radialgrid.setMajorRadius(R0)
    ss.radialgrid.setWallRadius(l_MK2)
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


def drawplot1(axs, so, color='r', toffset=0):
    """
    Draw a plot with a output from the given STREAMOutput object.
    """
    t = so.grid.t[:] + toffset
    plotInternal(axs[0], t[1:], np.diff(so.eqsys.W_cold[:,0]) / np.diff(so.grid.t[:]), ylabel=r'Power consumption (W/m$^3$)', color=color)
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
    plotInternal(axs[2], t[1:], P_net, ylbl, label='Net electron heating power', showlabel=showlabel, color='k', linestyle='-', ylim=[0, 2e5])

    axs[0].set_xlim([0, 0.05])
    axs[1].set_xlim([0, 0.05])
    axs[2].set_xlim([0, 0.05])
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
        ss11 = generate(prefill=5e-5, nt=10000)
        ss11.save(f'settings11{ext}.h5')
        so11 = runiface(ss11, f'output11{ext}.h5', quiet=False)

        ss12 = STREAMSettings(ss11)
        ss12.fromOutput(f'output11{ext}.h5')
        ss12.timestep.setTmax(0.1 - ss11.timestep.tmax)
        ss12.timestep.setNumberOfSaveSteps(0)
        ss12.timestep.setNt(1000)
        ss12.save(f'settings12{ext}.h5')
        so12 = runiface(ss12, f'output12{ext}.h5', quiet=False)
    else:
        so11 = STREAMOutput(f'output11{ext}.h5')
        so12 = STREAMOutput(f'output12{ext}.h5')

    if settings.skip is None or (len(settings.skip) > 0 and 2 not in settings.skip):
        print('RUN 2')
        ss21 = generate(prefill=7e-5, nt=10000)
        so21 = runiface(ss21, f'output21{ext}.h5', quiet=False)

        ss22 = STREAMSettings(ss21)
        ss22.fromOutput(f'output21{ext}.h5')
        ss22.timestep.setTmax(0.1 - ss21.timestep.tmax)
        ss22.timestep.setNumberOfSaveSteps(0)
        ss22.timestep.setNt(1000)
        so22 = runiface(ss22, f'output22{ext}.h5', quiet=False)
    else:
        so21 = STREAMOutput(f'output21{ext}.h5')
        so22 = STREAMOutput(f'output22{ext}.h5')

    if settings.plot:
        fig1, axs1 = plt.subplots(2, 1, figsize=(7,5), sharex=True)

        drawplot1(axs1, so11, color='b')
        drawplot1(axs1, so12, color='b', toffset=so11.grid.t[-1])
        drawplot1(axs1, so21, color='r')
        drawplot1(axs1, so22, color='r', toffset=so21.grid.t[-1])

        fig2, axs2 = plt.subplots(4, 1, figsize=(7, 10))
        drawplot2(axs2, so11, color='b')
        drawplot2(axs2, so12, color='b', toffset=so11.grid.t[-1])
        drawplot2(axs2, so21, color='r')
        drawplot2(axs2, so22, color='r', toffset=so21.grid.t[-1])

        axs2[0].legend([r'$p = 5\times 10^{-5}\,\mathrm{Torr}$', r'$p = 7\times 10^{-5}\,\mathrm{Torr}$'])

        fig3, axs3 = plt.subplots(3, 1, figsize=(7, 10))
        fig4, axs4 = plt.subplots(3, 1, figsize=(7, 10))

        drawplot3(axs3, so11)
        drawplot3(axs3, so12, toffset=so11.grid.t[-1], showlabel=False)

        drawplot3(axs4, so21)
        drawplot3(axs4, so22, toffset=so21.grid.t[-1], showlabel=False)

        plt.tight_layout()
        plt.show()

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


