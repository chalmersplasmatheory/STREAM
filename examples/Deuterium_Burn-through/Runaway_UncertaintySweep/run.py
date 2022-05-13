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



def drawplot4(axs, so, toffset=0, showlabel=True, save=False, first=True):
    """
    Draw a plot with a output from the given STREAMOutput object.
    """
    t = so.grid.t[:] + toffset


    Te = so.eqsys.T_cold[:, 0]
    Ti = so.eqsys.W_i.getTemperature()['D'][:, 0]

    ne = so.eqsys.n_cold[:, 0]
    nD0 = so.eqsys.n_i['D'][0][:]
    nD1 = so.eqsys.n_i['D'][1][:]

    vth = np.sqrt(2 * e * Te / m_e)
    streamingParameter = so.eqsys.j_tot[:, 0] / (e * ne * vth)

    Vp = so.other.stream.V_p[:, 0]
    Poh = -so.other.fluid.Tcold_ohmic[:, 0] * Vp
    Prad = so.other.fluid.Tcold_radiation[:, 0] * Vp
    Ptransp = so.other.stream.Tcold_transport[:, 0] * Vp
    Pequi = so.other.fluid.Tcold_ion_coll[:, 0] * Vp

    # Uext = so.eqsys.V
    Ures = 2 * np.pi * 5.7 * so.eqsys.E_field[:, 0]

    Ip = so.eqsys.I_p[:, 0]
    # Ire = e * c * 1.6**2 * np.pi * so.eqsys.n_re[:]
    Ire = so.eqsys.j_re.current()[:]
    Iohm = so.eqsys.j_ohm.current()[:]
    Iwall = so.eqsys.I_wall[:, 0]

    EoverED = so.eqsys.E_field.norm('ED')[1:, 0]
    ECoverED = so.other.fluid.Eceff[:, 0] / so.other.fluid.EDreic[:, 0]

    tauRE = so.other.stream.tau_RE[:, 0]
    tauRE1 = so.other.stream.tau_RE1[:, 0]
    tauRE2 = so.other.stream.tau_RE2[:, 0]

    totD = nD0[1:].flatten() * so.other.stream.V_n_tot['D'][:].flatten() + nD1[1:].flatten() * so.other.stream.V_p[:,
                                                                                               0].flatten()
    gammaD = (1 - nD0[1:].flatten() * so.other.stream.V_n_tot['D'][:].flatten() / totD)

    gammaTot = so.other.fluid.runawayRate[:, 0]
    gammaDreicer = so.other.fluid.gammaDreicer[:, 0]
    gammaAva = so.other.fluid.GammaAva[:, 0] * so.eqsys.n_re[1:, 0]

    if save:
        if first:
            worab = 'w'
        else:
            worab = 'ab'
        t_csv = open('Data/time.csv', worab)
        np.savetxt(t_csv, t)
        t_csv.close()
        Ip_csv = open('Data/PlasmaCurrent.csv', worab)
        np.savetxt(Ip_csv, Ip)
        Ip_csv.close()
        Ire_csv = open('Data/RunawayCurrent.csv', worab)
        np.savetxt(Ire_csv, Ire)
        Ire_csv.close()
        Iohm_csv = open('Data/OhmicCurrent.csv', worab)
        np.savetxt(Iohm_csv, Iohm)
        Iohm_csv.close()
        Te_csv = open('Data/ElectronTemperature.csv', worab)
        np.savetxt(Te_csv, Te)
        Te_csv.close()
        Ti_csv = open('Data/IonTemperature.csv', worab)
        np.savetxt(Ti_csv, Ti)
        Ti_csv.close()
        EoverED_csv = open('Data/ElectricField.csv', worab)
        np.savetxt( EoverED_csv,  EoverED)
        EoverED_csv.close()
        ECoverED_csv = open('Data/CriticalElectricField.csv', worab)
        np.savetxt(ECoverED_csv, ECoverED)
        ECoverED_csv.close()
        ne_csv = open('Data/ElectronDensity.csv', worab)
        np.savetxt(ne_csv, ne)
        ne_csv.close()

        gammaTot_csv = open('Data/TotalRunawayRate.csv', worab)
        np.savetxt(gammaTot_csv, gammaTot)
        gammaTot_csv.close()
        gammaDreicer_csv = open('Data/DreicerRunawayRate.csv', worab)
        np.savetxt(gammaDreicer_csv, gammaDreicer)
        gammaDreicer_csv.close()
        gammaAva_csv = open('Data/AvalancheRunawayRate.csv', worab)
        np.savetxt(gammaAva_csv, gammaAva)
        gammaAva_csv.close()
        tauRE_csv = open('Data/RunawayConfinementTime.csv', worab)
        np.savetxt(tauRE_csv, tauRE)
        tauRE_csv.close()
        tauRE1_csv = open('Data/RunawayParallellConfinementTime.csv', worab)
        np.savetxt(tauRE1_csv, tauRE1)
        tauRE1_csv.close()
        tauRE2_csv = open('Data/RunawayDriftConfinementTime.csv', worab)
        np.savetxt(tauRE2_csv, tauRE2)
        tauRE2_csv.close()

    plotInternal(axs[0, 0], t, Te, ylabel=r'$T$ (eV)', color='k', showlabel=showlabel, label=r'$T_{\rm e}$', yscalelog=True)
    plotInternal(axs[0, 0], t, Ti, ylabel=r'$T$ (eV)', color='m', showlabel=showlabel, label=r'$T_{\rm i}$', yscalelog=True)

    plotInternal(axs[0, 1], t, ne, ylabel=r'$n$ (m$^{-3}$)', color='k', showlabel=showlabel,
                 label=r'$n_{\rm e}$')
    plotInternal(axs[0, 1], t, nD0, ylabel=r'$n$ (m$^{-3}$)', color='r', showlabel=showlabel,
                 label=r'$n_{\rm D0}$')
    plotInternal(axs[0, 1], t, nD1, ylabel=r'$n$ (m$^{-3}$)', color='b', showlabel=showlabel,
                 label=r'$n_{\rm D1}$')

    plotInternal(axs[0, 2], t, streamingParameter, ylabel=r'$u_e/v_{\rm th}$', color='k', showlabel=showlabel,
                 label=r'$\xi$')
    axs[0, 2].set_ylim([-0.1 * np.max(streamingParameter[100:]), 1.1 * np.max(streamingParameter[100:])])

    plotInternal(axs[1, 0], t[1:], Poh / 1e6, ylabel=r'$P$ (MW)', color='m', showlabel=showlabel, label=r'$P_{\rm oh}$')
    plotInternal(axs[1, 0], t[1:], Prad / 1e6, ylabel=r'$P$ (MW)', color='r', showlabel=showlabel,
                 label=r'$P_{\rm rad}$')
    plotInternal(axs[1, 0], t[1:], Ptransp / 1e6, ylabel=r'$P$ (MW)', color='c', showlabel=showlabel,
                 label=r'$P_{\rm transp}$')
    plotInternal(axs[1, 0], t[1:], Pequi / 1e6, ylabel=r'$P$ (MW)', color='g', showlabel=showlabel,
                 label=r'$P_{\rm equi}$')

    # plotInternal(axs[1, 1], t, Uext, ylabel=r'$U$ (V)', color='k', showlabel=showlabel, label='$U_{\rm ext}$')
    plotInternal(axs[1, 1], t, Ures, ylabel=r'$U$ (V)', color='m', showlabel=showlabel, label=r'$U_{\rm res}$')

    plotInternal(axs[1, 2], t[1:], tauRE, ylabel=r'$\tau_{\rm re}$ (s)', color='k', showlabel=showlabel,
                 label=r'$\tau_{\rm re}$', yscalelog=True)
    plotInternal(axs[1, 2], t[1:], tauRE1, ylabel=r'$\tau_{\rm re}$ (s)', color='r', showlabel=showlabel,
                 label=r'Parallel', linestyle='--', yscalelog=True)
    plotInternal(axs[1, 2], t[1:], tauRE2, ylabel=r'$\tau_{\rm re}$ (s)', color='r', showlabel=showlabel,
                 label=r'Drifts', linestyle=':', yscalelog=True)

    plotInternal(axs[2, 0], t, Ip / 1e6, ylabel=r'$I$ (MA)', color='k', showlabel=showlabel, label=r'$I_{\rm p}$')
    plotInternal(axs[2, 0], t, Ire / 1e6, ylabel=r'$I$ (MA)', color='m', showlabel=showlabel, label=r'$I_{\rm re}$')
    plotInternal(axs[2, 0], t, Iohm / 1e6, ylabel=r'$I$ (MA)', color='r', showlabel=showlabel, label=r'$I_{\rm ohm}$')
    # axs[2, 0].set_ylim([0, 2])

    plotInternal(axs[2, 1], t[1:], EoverED * 100, ylabel=r'$E/E_{\mathrm{D}} (\%)$ ', color='k', showlabel=showlabel,
                 label=r'$E/E_{\rm D}$', yscalelog=False)
    plotInternal(axs[2, 1], t[1:], ECoverED * 100, ylabel=r'$E/E_{\mathrm{D}} (\%)$ ', color='m', showlabel=showlabel,
                 label=r'$E_{\rm C}/E_{\rm D}$', yscalelog=False)
    axs[2, 1].set_ylim([-0.1*np.max(EoverED[100:]*100), 1.1*np.max(EoverED[100:])*100])

    plotInternal(axs[2, 2], t, Iwall / 1e3, ylabel=r'$I_{\rm wall}$ (kA)', color='k', showlabel=showlabel,
                 label=r'$I_{\rm wall}$')

    plotInternal(axs[3, 0], t[1:], gammaD, ylabel=r'$\gamma$ (\%)', color='m', showlabel=showlabel,
                 label=r'$\gamma_{\rm D}$')

    plotInternal(axs[3, 1], t[1:], gammaTot, ylabel=r'$\gamma_{\rm re}$ (s$^{-1}$)', color='k', showlabel=showlabel,
                 label=r'$\mathrm{d}n_{\rm re} / \mathrm{d} t$', yscalelog=False)
    plotInternal(axs[3, 1], t[1:], gammaDreicer, ylabel=r'$\gamma_{\rm re}$ (s$^{-1}$)', color='b', showlabel=showlabel,
                 label=r'$\gamma_{\rm Dreicer}$', yscalelog=False)
    plotInternal(axs[3, 1], t[1:], gammaAva, ylabel=r'$\gamma_{\rm re}$ (s$^{-1}$)', color='r', showlabel=showlabel,
                 label=r'$\gamma_{\rm ava}$', yscalelog=False)
    #axs[3, 1].set_ylim([-0.1*np.max(gammaTot[10:]), 1.1*np.max(gammaTot[10:])])

    for i in range(axs.shape[0]):
        for j in range(axs.shape[1]):
            axs[i, j].set_xlim([0, t[-1]])
            axs[i, j].grid(True)


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

def makeplots(so11, so12, save=False):
    '''
    fig1, axs1 = plt.subplots(2, 1, figsize=(7,5), sharex=True)
    drawplot1(axs1, so11, color='b')
    drawplot1(axs1, so12, color='b', toffset=so11.grid.t[-1])
    fig2, axs2 = plt.subplots(4, 1, figsize=(7, 10))
    drawplot2(axs2, so11, color='b')
    drawplot2(axs2, so12, color='b', toffset=so11.grid.t[-1])
    fig3, axs3 = plt.subplots(3, 1, figsize=(7, 10))
    drawplot3(axs3, so11)
    drawplot3(axs3, so12, toffset=so11.grid.t[-1], showlabel=False)
    '''
    fig4, axs4 = plt.subplots(4, 3, figsize=(12, 10))
    drawplot4(axs4, so11, showlabel=True, save=save)
    drawplot4(axs4, so12, toffset=so11.grid.t[-1], showlabel=False, save=save, first=False)

    axs4[0, 0].legend(frameon=False, prop={'size': 10})
    axs4[0, 1].legend(frameon=False, prop={'size': 10})
    axs4[1, 0].legend(frameon=False, prop={'size': 10})
    axs4[2, 0].legend(frameon=False, prop={'size': 10})
    axs4[2, 1].legend(frameon=False, prop={'size': 10})
    axs4[1, 2].legend(frameon=False, prop={'size': 10})
    axs4[3, 1].legend(frameon=False, prop={'size': 10})

    plt.tight_layout()
    plt.show()

def saveDataArticle():
    pgp = 2 / 133.32 * 8e-5
    nPoints = 100
    uncertaintyFactors = np.logspace(-3, 0, nPoints)
    nre_peak = np.zeros(nPoints)
    t_thresh = np.zeros(nPoints)
    loss = np.zeros(nPoints)
    loss_ava = np.zeros(nPoints)
    for f_u, i in zip(uncertaintyFactors, range(nPoints)):
        ss21 = generate(prefill=pgp, uncertaintyFactor=f_u)
        ss21.save(f'SweepSettings/settings1_{f_u}.h5')
        so21 = runiface(ss21, f'SweepSettings/output1_{f_u}.h5', quiet=False)

        ss22 = STREAMSettings(ss21)
        ss22.fromOutput(f'SweepSettings/output1_{f_u}.h5')
        ss22.timestep.setTmax(0.15 - ss21.timestep.tmax)
        ss22.timestep.setNt(1000)
        ss22.save(f'SweepSettings/settings2_{f_u}.h5')
        so22 = runiface(ss22, f'SweepSettings/output2_{f_u}.h5', quiet=False)

        nre_peak[i], t_thresh[i], loss[i], loss_ava[i] = getNrePeak(so22, so21, f_u, 9e13)
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
