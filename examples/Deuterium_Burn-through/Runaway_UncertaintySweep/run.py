#!/usr/bin/env python3
#
# This example uses parameters from Table 4.1 of [H. T. Kim, PhD Thesis (2013)]
# to simulate basic Deuterium burn-through in a JET-like plasma.
#

import argparse
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.constants import c, e, m_e, mu_0
from scipy.signal import find_peaks
sys.path.append('../../../py')

import PlasmaParameters as Formulas
import DREAM.Settings.Equations.ColdElectronTemperature as Tcold
import DREAM.Settings.Solver as Solver
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Equations.ElectricField as ElectricField_D
import DREAM.Settings.Atomics as Atomics
from STREAM import STREAMOutput, STREAMSettings, runiface
import STREAM.Settings.Equations.ElectricField as ElectricField
import STREAM.Settings.Equations.IonSpecies as Ions


def generate(prefill=1e-6, gamma=2e-3, Vloop=12, Vloop_t=0, I0=2.4e3, tmax=1e-4, nt=4000, EfieldDYON=True, uncertaintyFactor=1.0):
    """
    Generate a STREAMSettings object for a simulation with the specified
    parameters.
    
    :param prefill:   Prefill gas pressure [Torr]
    :param gamma:     Degree of ionization
    :param fractionT: Fraction of tritium
    :param Vloop:     External loop voltage [V]
    :param Vloop_t:   Time vector corresponding to given loop voltage [V]
    :param I0:        Initial plasma current [A]
    :param tmax:      Simulation time [s]
    :param nt:        Number of time steps
    """

    n0 = 3.22e22 * prefill  # Initial total deuterium density
    nD = n0 * np.array([[1-gamma], [gamma]])

    a = 1.6  # Plasma minor radius [m]
    R0 = 5.65  # Plasma major radius [m]
    Btor = 2.65 * 6.2 / R0  # Toroidal magnetic field [T]
    l_MK2 = 1       # Distance between plasma centre and passive structure [m]
    V_vessel = 1000  # Vacuum vessel volume

    Te0 = 1     # electron temperature [eV]
    Ti0 = 0.026 # ion temperature [eV]

    # Initial electric field
    j0 = I0 / (a ** 2 * np.pi)
    E0 = j0 / Formulas.evaluateSpitzerConductivity(n=nD[1], T=Te0, Z=1)

    ss = STREAMSettings()

    ss.atomic.adas_interpolation = Atomics.ADAS_INTERP_BILINEAR

    # Electric field
    if EfieldDYON:
        Lp = float(mu_0 * R0 * (np.log(8 * R0 / a) + 0.25 - 2))

        ss.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
        ss.eqsys.E_field.setInitialProfile(E0)
        ss.eqsys.E_field.setInductances(Lp=Lp, Lwall=9.1e-6, M=2.49e-6, Rwall=1e7)
        ss.eqsys.E_field.setCircuitVloop(Vloop, Vloop_t)
    else:
        R = 1e7
        L = 9.1e-6
        wall_time = L / R

        ss.eqsys.E_field.setType(ElectricField_D.TYPE_SELFCONSISTENT)
        ss.eqsys.E_field.setInitialProfile(efield=E0)
        ss.eqsys.E_field.setBoundaryCondition(ElectricField_D.BC_TYPE_TRANSFORMER, V_loop_wall_R0=Vloop / R0,
                                              times=Vloop_t, inverse_wall_time=1 / wall_time, R0=R0)

    # Electron temperature
    ss.eqsys.T_cold.setType(Tcold.TYPE_SELFCONSISTENT)
    ss.eqsys.T_cold.setInitialProfile(Te0)
    ss.eqsys.T_cold.setRecombinationRadiation(recombination=Tcold.RECOMBINATION_RADIATION_INCLUDED)

    # Ions
    ss.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=nD, r=np.array([0]), T=Ti0)

    # Enable runaway
    ss.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
    ss.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)

    # Recycling coefficients (unused)
    ss.eqsys.n_i.setJET_CWrecycling()

    # Radial grid
    ss.radialgrid.setB0(Btor)
    ss.radialgrid.setMinorRadius(a)
    ss.radialgrid.setMajorRadius(R0)
    ss.radialgrid.setWallRadius(l_MK2)
    ss.radialgrid.setVesselVolume(V_vessel)
    ss.radialgrid.setRunawayConfinementUncertainty(uncertaintyFactor)

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

    ss.other.include('fluid', 'stream', 'scalar')

    return ss

def getNrePeak(so, so_init, factor, thresh, toffset=0):
    t1 = so_init.grid.t[:]
    dt1 = (t1[-1] - t1[0]) / (len(t1) - 1)
    dreicer1 = so_init.other.fluid.gammaDreicer[:, 0]
    dreicer_int1 = np.sum(dreicer1 * dt1)
    t2 = so.grid.t[:] + 1e-4
    dt2 = (t2[-1] - t2[0]) / (len(t2) - 1)
    dreicer2 = so.other.fluid.gammaDreicer[:, 0]
    dreicer_int2 = np.sum(dreicer2 * dt2)
    loss = 1-so.eqsys.n_re[-1][0]/(dreicer_int1+dreicer_int2)
    ava1 = so_init.other.fluid.GammaAva[:, 0] * so_init.eqsys.n_re[1:, 0]
    ava2 = so.other.fluid.GammaAva[:, 0] * so.eqsys.n_re[1:, 0]
    ava_int = np.sum(ava1*dt1) + np.sum(ava2*dt2)
    loss_ava = 1 - so.eqsys.n_re[-1][0] / (dreicer_int1 + dreicer_int2 + ava_int)

    n_re = so.eqsys.n_re[:].flatten()
    t = so.grid.t[:] + toffset
    iPeaks = find_peaks(n_re)[0]

    iThresh = np.abs(n_re[iPeaks[0]:]-thresh).argmin() + iPeaks[0]
    if iThresh == len(n_re)-1 or n_re[iThresh+1] > thresh:
        t[iThresh] = np.inf

    #if len(iPeaks)!=1:
    #plt.plot(t,n_re)
    #plt.plot(np.array([t[0], t[-1]]), np.array([thresh, thresh]))
    #plt.scatter(t[iPeaks], n_re[iPeaks])
    #if iThresh != -1:
    #plt.scatter(t[iThresh], n_re[iThresh])
    #plt.title('Uncertainty factor is ' + str(factor))
    #plt.show()
    return n_re[iPeaks[0]], t[iThresh] - t[iPeaks[0]], loss, loss_ava



def drawplot4(axs, so, toffset=0, showlabel=True, save=False, first=True, fileaddon='', color='k', linestyle='solid', label=''):
    """
    Draw a plot with a output from the given STREAMOutput object.
    """
    t = so.grid.t[:] + toffset

    nre = so.eqsys.n_re[:].flatten()

    if save:
        if first:
            worab = 'w'
        else:
            worab = 'ab'
        t_csv = open('Data/time' + fileaddon + '.csv', worab)
        np.savetxt(t_csv, t)
        t_csv.close()
        nre_csv = open('Data/RunawayElectronDensity' + fileaddon + '.csv', worab)
        np.savetxt(nre_csv, nre)
        nre_csv.close()


    plotInternal(axs, t, nre, ylabel=r'$n_{\rm re}$ (m$^{-3}$)', color=color, linestyle=linestyle, showlabel=showlabel, label=label)

    axs.grid(True, color='gainsboro', linewidth=0.5)


def plotInternal(ax, x, y, ylabel, xlbl=True, ylim=None, log=False, showlabel=False, label=None, yscalelog = False, *args, **kwargs):
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

    if yscalelog:
        ax.set_yscale('log')

def makeplot(axs, so1, so2, save=False, fileaddon='', color='k', linestyle='solid', label=''):
    drawplot4(axs, so1, showlabel=True, save=save, fileaddon=fileaddon, color=color, linestyle=linestyle, label=label)
    drawplot4(axs, so2, toffset=so1.grid.t[-1], showlabel=False, save=save, first=False, fileaddon=fileaddon, color=color, linestyle=linestyle, label=label)



def saveDataArticle():
    pgp = 2 / 133.32 * 8e-5
    nPoints = 100
    uncertaintyFactors = [1, 0.1, 1e-3]#np.logspace(-3, 0, nPoints)
    nre_peak = np.zeros(nPoints)
    t_thresh = np.zeros(nPoints)
    loss = np.zeros(nPoints)
    loss_ava = np.zeros(nPoints)

    fig, axs = plt.subplots(1, 1, figsize=(5, 5))
    fileaddons = ['_1', '_0.1', '_0.001']
    colors = ['tab:blue', 'k', 'tab:red']
    linestyles = ['solid', 'dashed', 'dotted']
    labels = ['$f_u=1$', '$f_u=0.1$', '$f_u=0.001$']
    for f_u, i, fa, c, ls, lab in zip(uncertaintyFactors, range(nPoints), fileaddons, colors, linestyles, labels):
        print(str(f_u))
        ss1 = generate(prefill=pgp, uncertaintyFactor=f_u)
        ss1.save(f'SweepSettings/settings1_{f_u}.h5')
        so1 = runiface(ss1, f'SweepSettings/output1_{f_u}.h5', quiet=False)

        ss2 = STREAMSettings(ss1)
        ss2.fromOutput(f'SweepSettings/output1_{f_u}.h5')
        ss2.timestep.setTmax(0.15 - ss1.timestep.tmax)
        ss2.timestep.setNt(1000)
        ss2.save(f'settings2_{f_u}.h5')
        so2 = runiface(ss2, f'output2_{f_u}.h5', quiet=False)

        makeplot(axs, so1, so2, save=True, fileaddon=fa, color=c, linestyle=ls, label=lab)
    plt.legend()
    plt.tight_layout()
    plt.show()

        #nre_peak[i], t_thresh[i], loss[i], loss_ava[i] = getNrePeak(so2, so1, f_u, 9e13)
    '''
    plt.semilogx(uncertaintyFactors, nre_peak)
    plt.show()
    plt.semilogx(uncertaintyFactors, t_thresh)
    plt.show()
    plt.semilogx(uncertaintyFactors, loss)
    plt.show()
    plt.semilogx(uncertaintyFactors, loss_ava)
    plt.show()
    fu_csv = open('Data/uncertaintyfactor.csv', "w")
    np.savetxt(fu_csv, uncertaintyFactors)
    fu_csv.close()
    nrep_csv = open('Data/peakvalueNre.csv', "w")
    np.savetxt(nrep_csv, nre_peak)
    nrep_csv.close()
    tt_csv = open('Data/tThreshold.csv', "w")
    np.savetxt(tt_csv, t_thresh)
    tt_csv.close()
    loss_csv = open('Data/RunawayLoss.csv', "w")
    np.savetxt(loss_csv, loss)
    loss_csv.close()
    lossava_csv = open('Data/RunawayLossWithAvalanche.csv', "w")
    np.savetxt(lossava_csv, loss_ava)
    lossava_csv.close()
    '''

def main(argv):
    FONTSIZE = 16
    plt.rcParams.update({'font.size': FONTSIZE})
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--skip-all', help="Skip all simulations and only generate plots", dest="skip", action='store_const', const=[1,2])
    parser.add_argument('-e', '--extension', help="Append the specified string to the end of file names", dest="extension", action='store', default='')
    parser.add_argument('-n', '--no-plot', help="Only run simulations, don't generate plots", dest="plot", action='store_false', default=True)
    parser.add_argument('-s', '--skip', help="Skip the specified simulation(s)", dest="skip", nargs="*", type=int)
    settings = parser.parse_args()
    ext = '' if not settings.extension else '_' + settings.extension
    if settings.skip is None or (len(settings.skip) > 0 and 2 not in settings.skip):
        ss21 = generate()
        ss21.save(f'settings1WithT{ext}.h5')
        so21 = runiface(ss21, f'output1WithT{ext}.h5', quiet=False)
        ss22 = STREAMSettings(ss21)
        ss22.fromOutput(f'output1WithT{ext}.h5')
        ss22.timestep.setTmax(8 - ss21.timestep.tmax)
        ss22.timestep.setNt(20000)
        ss22.save(f'settings2WithT{ext}.h5')
        so22 = runiface(ss22, f'output2WithT{ext}.h5', quiet=False)
    else:
        so21 = STREAMOutput(f'output1WithT{ext}.h5')
        so22 = STREAMOutput(f'output2WithT{ext}.h5')
    if settings.plot:
        makeplots(so21, so22)
    '''
    saveDataArticle()

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
