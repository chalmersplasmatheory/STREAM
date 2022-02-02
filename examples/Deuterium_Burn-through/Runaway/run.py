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


def generate(prefill=1e-6, gamma=2e-2, Vloop=12, Vloop_t=0, I0=40e3, tmax=1e-4, nt=10000, EfieldDYON=True, runaways=True):
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

    Btor = 2.65     # Toroidal magnetic field [T]
    a = 1.6         # Plasma minor radius [m]
    R0 = 5.65       # Plasma major radius [m]
    l_MK2 = 1       # Distance between plasma centre and passive structure [m]
    V_vessel = 1700  # Vacuum vessel volume

    Te0 = 5     # electron temperature [eV]
    Ti0 = 1 # ion temperature [eV]

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
    if runaways:
        ss.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
        ss.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)
    else:
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

def drawplot4(axs, so, toffset=0, showlabel=True, save=False, fileaddon=''):
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
        t_csv = open('Data/time_'+ fileaddon +'.csv', 'ab')
        np.savetxt(t_csv, t)
        t_csv.close()
        Ip_csv = open('Data/PlasmaCurrent_'+ fileaddon +'.csv', 'ab')
        np.savetxt(Ip_csv, Ip)
        Ip_csv.close()
        Ire_csv = open('Data/RunawayCurrent_' + fileaddon + '.csv', 'ab')
        np.savetxt(Ire_csv, Ire)
        Ire_csv.close()
        Iohm_csv = open('Data/OhmicCurrent_' + fileaddon + '.csv', 'ab')
        np.savetxt(Iohm_csv, Iohm)
        Iohm_csv.close()
        Te_csv = open('Data/ElectronTemperature_'+ fileaddon +'.csv', 'ab')
        np.savetxt(Te_csv, Te)
        Te_csv.close()
        Ti_csv = open('Data/IonTemperature_'+ fileaddon +'.csv', 'ab')
        np.savetxt(Ti_csv, Ti)
        Ti_csv.close()
        EoverED_csv = open('Data/ElectricField_'+ fileaddon +'.csv', 'ab')
        np.savetxt( EoverED_csv,  EoverED)
        EoverED_csv.close()
        ECoverED_csv = open('Data/CriticalElectricField_' + fileaddon + '.csv', 'ab')
        np.savetxt(ECoverED_csv, ECoverED)
        ECoverED_csv.close()

        gammaTot_csv = open('Data/TotalRunawayRate_' + fileaddon + '.csv', 'ab')
        np.savetxt(gammaTot_csv, gammaTot)
        gammaTot_csv.close()
        gammaDreicer_csv = open('Data/DreicerRunawayRate_' + fileaddon + '.csv', 'ab')
        np.savetxt(gammaDreicer_csv, gammaDreicer)
        gammaDreicer_csv.close()
        gammaAva_csv = open('Data/AvalancheRunawayRate_' + fileaddon + '.csv', 'ab')
        np.savetxt(gammaAva_csv, gammaAva)
        gammaAva_csv.close()
        tauRE_csv = open('Data/RunawayConfinementTime_'+ fileaddon +'.csv', 'ab')
        np.savetxt(tauRE_csv, tauRE)
        tauRE_csv.close()
        tauRE1_csv = open('Data/RunawayDriftConfinementTime_' + fileaddon + '.csv', 'ab')
        np.savetxt(tauRE1_csv, tauRE1)
        tauRE1_csv.close()
        tauRE2_csv = open('Data/RunawayParallellConfinementTime_' + fileaddon + '.csv', 'ab')
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

def makeplots(so11, so12, save=False, fileaddon=''):
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
    drawplot4(axs4, so11, showlabel=True, save=save, fileaddon=fileaddon)
    drawplot4(axs4, so12, toffset=so11.grid.t[-1], showlabel=False, save=save, fileaddon=fileaddon)

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
    fileaddons = ['noRE', 'enabledRE', 'disabledRE']
    runaways = [True, True, False]
    prefillpressures = 2 / 133.32 * np.array([7e-4, 7e-5, 7e-5]) # Pa -> Torr
    for fa, re, pgp in zip(fileaddons, runaways, prefillpressures):
        ss21 = generate(prefill=pgp, runaways=re)
        ss21.save(f'settings1WithT{fa}.h5')
        so21 = runiface(ss21, f'output1WithT{fa}.h5', quiet=False)

        ss22 = STREAMSettings(ss21)
        ss22.fromOutput(f'output1WithT{fa}.h5')
        ss22.timestep.setTmax(8 - ss21.timestep.tmax)
        ss22.timestep.setNt(100000)
        ss22.save(f'settings2WithT{fa}.h5')
        so22 = runiface(ss22, f'output2WithT{fa}.h5', quiet=False)

        makeplots(so21, so22, save=True, fileaddon=fa)

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


