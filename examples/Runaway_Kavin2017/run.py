#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, e, m_e, mu_0
import sys
sys.path.append('../Deuterium_Burn-through')
sys.path.append('../../py')

#from run import makeplots
import DREAM.Settings.Equations.ColdElectronTemperature as Tcold
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver

import DREAM.Settings.Atomics as Atomics
from STREAM import STREAMOutput, STREAMSettings, runiface
import STREAM.Settings.Equations.ElectricField as ElectricField
import DREAM.Settings.Equations.ElectricField as ElectricField_D
import STREAM.Settings.Equations.IonSpecies as Ions



def generate(n0=0.01e20/7, Btor = 2.65, gamma=3e-2, Vloop=10.6, Vloop_t=0, Ures0=14, tmax=0.003, nt=10000, Rwall=1e7, EfieldDyon=False, tritium=False, impurity=True):
    """
    Generate a STREAMSettings object for a simulation with the specified
    parameters.

    floats:
    :param n0:      Deuterium density [m^-3]
    :param Btor:    Toroidal magnetic field [T]
    :param gamma:   Degree of ionization
    :param Vloop:   External loop voltage [V]
    :param Vloop_t: Time vector corresponding to given loop voltage [V]
    :param Ures0:   Initial voltage in plasma [V]
    :param tmax:    Simulation time [s]
    :param nt:      Number of  timesteps
    :param Rwall:   Wall resistance [ohm]
    :param nt:      Number of time steps

    booleans:
    :param EfieldDyon:  If circuit equations from DYON should be used, otherwise selfconsistent evolution of electric field
    :param tritium:     If half of the hydrogen plasma should be tritium, otherwise only deuterium
    :param impurity:    If there should be impurities, in this case iron
    """
    if tritium:
        nD = np.array([[n0*0.5], [n0 * gamma / (1 - gamma)*0.5]])
        nT = np.array([[n0 * 0.5], [n0 * gamma / (1 - gamma) * 0.5]])
    else:
        nD = np.array([[n0], [n0 * gamma / (1 - gamma)]])

    a = 1.6         # Plasma minor radius [m]
    R0 = 5.7        # Plasma major radius [m]
    b = 8*R0/np.exp(7/4)
    V_vessel = 1700 # Vacuum vessel volume [m^3]

    Te0 = 5     # electron temperature [eV]
    Ti0 = 1     # ion temperature [eV]

    impurity = False


    # Initial electric field
    E0 = Ures0/(2 * np.pi * R0)

    ss = STREAMSettings()

    ss.atomic.adas_interpolation = Atomics.ADAS_INTERP_BILINEAR

    # Electric field
    Lwall = 9.1e-6
    if EfieldDyon:
        ss.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
        ss.eqsys.E_field.setInitialProfile(E0)
        Lp = float(mu_0 * R0 * (np.log(8 * R0 / a) + 0.25 - 2))
        #ss.eqsys.E_field.setInductances(Lp=Lp, Lwall=9.1e-6, M=2.49e-6, Rwall=5e-6)
        ss.eqsys.E_field.setInductances(Lp=Lp, Lwall=Lwall, M=2.49e-6, Rwall=Rwall)
        ss.eqsys.E_field.setCircuitVloop(Vloop, Vloop_t)
        # Value based on statement in (de Vries & Gribov 2019; page 12)
        #ss.eqsys.E_field.setWallCurrent(1.5e6)
    else:
        ss.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
        ss.eqsys.E_field.setInitialProfile(efield=E0)

        walltime = Lwall / Rwall
        ss.eqsys.E_field.setBoundaryCondition(ElectricField_D.BC_TYPE_TRANSFORMER, V_loop_wall_R0=Vloop/R0, times=Vloop_t, inverse_wall_time=1/walltime, R0=R0)

    # Electron temperature
    ss.eqsys.T_cold.setType(Tcold.TYPE_SELFCONSISTENT)
    ss.eqsys.T_cold.setInitialProfile(Te0)
    ss.eqsys.T_cold.setRecombinationRadiation(recombination=Tcold.RECOMBINATION_RADIATION_INCLUDED)

    # Ions
    ss.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=nD, r=np.array([0]), T=Ti0)
    if tritium:
        ss.eqsys.n_i.addIon(name='T', Z=1, iontype=Ions.IONS_DYNAMIC, n=nT, r=np.array([0]), T=Ti0)

    if impurity:
        ss.eqsys.n_i.addIon(name='Fe', Z=26, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=1e3, r=np.array([0]), T=Ti0)

    # Enable runaway
    ss.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
    ss.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)

    # Recycling coefficients
    iD = ss.eqsys.n_i.getIndex('D')
    ss.eqsys.n_i.ions[iD].setRecyclingCoefficient('D', 1)
    if impurity:
        ss.eqsys.n_i.ions[iD].setRecyclingCoefficient('Fe', 0.00001)
    if tritium:
        iT = ss.eqsys.n_i.getIndex('T')
        ss.eqsys.n_i.ions[iT].setRecyclingCoefficient('T', 1)
        if impurity:
            ss.eqsys.n_i.ions[iT].setRecyclingCoefficient('Fe', 0.00001)


    # Radial grid
    ss.radialgrid.setB0(Btor)
    ss.radialgrid.setMinorRadius(a)
    ss.radialgrid.setMajorRadius(R0)
    ss.radialgrid.setWallRadius(b)
    ss.radialgrid.setVesselVolume(V_vessel)
    ss.radialgrid.setIref(1e5)
    ss.radialgrid.setBv(2.0e-3)
    ss.radialgrid.setConnectionLengthFactor(3.0)

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


def drawplot1(axs, so, toffset=0.7, showlabel=True):
    """
    Draw a plot with a output from the given STREAMOutput object.
    """
    t = so.grid.t[:] + toffset

    Te = so.eqsys.T_cold[:, 0]
    Ti = so.eqsys.W_i.getTemperature()['D'][:, 0]

    ne = so.eqsys.n_cold[:, 0]
    nD0 = so.eqsys.n_i['D'][0][:]
    nD1 = so.eqsys.n_i['D'][1][:]

    vth = np.sqrt(2*e*Te / m_e)
    streamingParameter = so.eqsys.j_tot[:,0] / (e * ne * vth)

    nFe = []
    for iFe in range(0,27):
        nFe.append(so.eqsys.n_i['Fe'][iFe][:])
    nFe = np.array(nFe)

    totFe = nFe[0,1:].flatten() * so.other.stream.V_n_tot['Fe'][:].flatten()
    for iFe in range(0,27):
        totFe = totFe + nFe[iFe,1:].flatten() * so.other.stream.V_p[:, 0].flatten()

    gammaFe = (1 - nFe[0,1:].flatten() * so.other.stream.V_n_tot['Fe'][:].flatten() / totFe)

    totD = nD0[1:].flatten() * so.other.stream.V_n_tot['D'][:].flatten() + nD1[1:].flatten() * so.other.stream.V_p[:, 0].flatten()
    gammaD = (1 - nD0[1:].flatten() * so.other.stream.V_n_tot['D'][:].flatten() / totD)

    Vp = so.other.stream.V_p[:, 0]
    Poh = -so.other.fluid.Tcold_ohmic[:, 0] * Vp
    Prad = so.other.fluid.Tcold_radiation[:, 0] * Vp
    Ptransp = so.other.stream.Tcold_transport[:,0] * Vp
    Pequi = so.other.fluid.Tcold_ion_coll[:, 0] * Vp

    Ures = 2 * np.pi * 5.7 * so.eqsys.E_field[:,0]

    Ip = so.eqsys.I_p[:,0]
    Ire = so.eqsys.j_re.current()[:]
    Iwall = so.eqsys.I_wall[:,0]

    EoverED = so.eqsys.E_field.norm('ED')[1:,0]
    ECoverED = so.other.fluid.Eceff[:,0] / so.other.fluid.EDreic[:,0]

    tauRE = so.other.stream.tau_RE[:,0]
    tauRE1 = so.other.stream.tau_RE1[:,0]
    tauRE2 = so.other.stream.tau_RE2[:,0]

    gammaTot = so.other.fluid.runawayRate[:,0]
    gammaDreicer = so.other.fluid.gammaDreicer[:,0]
    gammaAva = so.other.fluid.GammaAva[:,0] * so.eqsys.n_re[1:,0]

    plotInternal(axs[0, 0], t, Te / 1e3, ylabel=r'$T$ (keV)', color='k', showlabel=showlabel, label=r'$T_{\rm e}$')
    plotInternal(axs[0, 0], t, Ti / 1e3, ylabel=r'$T$ (keV)', color='m', showlabel=showlabel, label=r'$T_{\rm i}$')

    plotInternal(axs[0, 1], t, ne / 1e20, ylabel=r'$n$ (1e20 m$^{-3}$)', color='k', showlabel=showlabel, label=r'$n_{\rm e}$')
    plotInternal(axs[0, 1], t, 7 * nD0 / 1e20, ylabel=r'$n$ (1e20 m$^{-3}$)', color='m', showlabel=showlabel, label=r'$n_{\rm D0}$')
    plotInternal(axs[0, 1], t, nD1 / 1e20, ylabel=r'$n$ (1e20 m$^{-3}$)', color='c', showlabel=showlabel, label=r'$n_{\rm D1}$')

    plotInternal(axs[0, 2], t, streamingParameter, ylabel=r'$u_e/v_{\rm th}$', color='k', showlabel=showlabel, label=r'$\xi$')
    axs[0,2].set_ylim([0, 0.2])

    plotInternal(axs[1, 0], t[1:], Poh / 1e6, ylabel=r'$P$ (MW)', color='m', showlabel=showlabel, label=r'$P_{\rm oh}$')
    plotInternal(axs[1, 0], t[1:], Prad / 1e6, ylabel=r'$P$ (MW)', color='r', showlabel=showlabel, label=r'$P_{\rm rad}$')
    plotInternal(axs[1, 0], t[1:], Ptransp / 1e6, ylabel=r'$P$ (MW)', color='c', showlabel=showlabel, label=r'$P_{\rm transp}$')
    plotInternal(axs[1, 0], t[1:], Pequi / 1e6, ylabel=r'$P$ (MW)', color='g', showlabel=showlabel, label=r'$P_{\rm equi}$')

    #plotInternal(axs[1, 1], t, Uext, ylabel=r'$U$ (V)', color='k', showlabel=showlabel, label='$U_{\rm ext}$')
    plotInternal(axs[1, 1], t, Ures, ylabel=r'$U$ (V)', color='m', showlabel=showlabel, label=r'$U_{\rm res}$')

    plotInternal(axs[1, 2], t[1:], tauRE, ylabel=r'$\tau_{\rm re}$ (s)', color='k', showlabel=showlabel, label=r'$\tau_{\rm re}$', yscalelog=True)
    plotInternal(axs[1, 2], t[1:], tauRE1, ylabel=r'$\tau_{\rm re}$ (s)', color='r', showlabel=showlabel, label=r'Parallel', linestyle='--', yscalelog=True)
    plotInternal(axs[1, 2], t[1:], tauRE2, ylabel=r'$\tau_{\rm re}$ (s)', color='r', showlabel=showlabel, label=r'Drifts', linestyle=':', yscalelog=True)

    plotInternal(axs[2, 0], t, Ip / 1e6, ylabel=r'$I$ (MA)', color='k', showlabel=showlabel, label=r'$I_{\rm p}$')
    plotInternal(axs[2, 0], t, Ire / 1e6, ylabel=r'$I$ (MA)', color='m', showlabel=showlabel, label=r'$I_{\rm re}$')
    #axs[2, 0].set_ylim([0, 2])

    plotInternal(axs[2, 1], t[1:], EoverED*100, ylabel=r'$E/E_{\mathrm{D}}$ (\%)', color='k', showlabel=showlabel, label=r'$E/E_{\rm D}$', yscalelog = False)
    plotInternal(axs[2, 1], t[1:], ECoverED*100, ylabel=r'$E/E_{\mathrm{D}}$ (\%)', color='m', showlabel=showlabel, label=r'$E_{\rm C}/E_{\rm D}$', yscalelog = False)
    axs[2,1].set_ylim([0, 10])

    plotInternal(axs[2, 2], t, Iwall/1e3, ylabel=r'$I_{\rm wall}$ (kA)', color='k', showlabel=showlabel, label=r'$I_{\rm wall}$')

    plotInternal(axs[3, 0], t[1:], gammaD, ylabel=r'$\gamma_{\rm D}$ (\%)', color='m', showlabel=showlabel, label=r'$\gamma_{\rm D}$')
    plotInternal(axs[3, 0], t[1:], gammaFe, ylabel=r'$\gamma_{\rm Fe}$ (\%)', color='k', showlabel=showlabel, label=r'$\gamma_{\rm Fe}$')
    
    n_norm = 1e20
    plotInternal(axs[3, 1], t[:], nFe[3, :]/n_norm, ylabel=r'$n$ m$^{-3}$', color='k', showlabel=showlabel, label=r'$n_{\rm Fe3}$', yscalelog = True)
    plotInternal(axs[3, 1], t[:], nFe[4, :]/n_norm, ylabel=r'$n$ m$^{-3}$', color='y', showlabel=showlabel, label=r'$n_{\rm Fe4}$', yscalelog = True)
    plotInternal(axs[3, 1], t[:], nFe[5, :]/n_norm, ylabel=r'$n$ m$^{-3}$', color='m', showlabel=showlabel, label=r'$n_{\rm Fe5}$', yscalelog = True)
    plotInternal(axs[3, 1], t[:], nFe[6, :]/n_norm, ylabel=r'$n$ m$^{-3}$', color='r', showlabel=showlabel, label=r'$n_{\rm Fe6}$', yscalelog = True)
    plotInternal(axs[3, 1], t[:], nFe[7, :]/n_norm, ylabel=r'$n$ m$^{-3}$', color='c', showlabel=showlabel, label=r'$n_{\rm Fe7}$', yscalelog = True)
    plotInternal(axs[3, 1], t[:], nFe[8, :]/n_norm, ylabel=r'$n$ m$^{-3}$', color='g', showlabel=showlabel, label=r'$n_{\rm Fe8}$', yscalelog = True)
    plotInternal(axs[3, 1], t[:], nFe[9, :]/n_norm, ylabel=r'$n$ m$^{-3}$', color='b', showlabel=showlabel, label=r'$n_{\rm Fe9}$', yscalelog = True)
    plotInternal(axs[3, 1], t[:], nFe[10, :]/n_norm, ylabel=r'$n$ m$^{-3}$', color='silver', showlabel=showlabel, label=r'$n_{\rm Fe10}$', yscalelog = True)
    plotInternal(axs[3, 1], t[:], nFe[11, :]/n_norm, ylabel=r'$n$ m$^{-3}$', color='dimgrey', showlabel=showlabel, label=r'$n_{\rm Fe11}$', yscalelog = True)
    plotInternal(axs[3, 1], t[:], nFe[12, :]/n_norm, ylabel=r'$n$ m$^{-3}$', color='gold', showlabel=showlabel, label=r'$n_{\rm Fe12}$', yscalelog = True)
    plotInternal(axs[3, 1], t[:], nFe[13, :]/n_norm, ylabel=r'$n$ m$^{-3}$', color='deeppink', showlabel=showlabel, label=r'$n_{\rm Fe13}$', yscalelog = True)

    plotInternal(axs[3, 2], t[1:], gammaTot, ylabel=r'$\gamma_{\rm re}$ (s$^{-1}$)', color='k', showlabel=showlabel, label=r'$\mathrm{d}n_{\rm re} / \mathrm{d} t$', yscalelog=False)
    plotInternal(axs[3, 2], t[1:], gammaDreicer, ylabel=r'$\gamma_{\rm re}$ (s$^{-1}$)', color='b', showlabel=showlabel, label=r'$\gamma_{\rm Dreicer}$', yscalelog=False)
    plotInternal(axs[3, 2], t[1:], gammaAva, ylabel=r'$\gamma_{\rm re}$ (s$^{-1}$)', color='r', showlabel=showlabel, label=r'$\gamma_{\rm ava}$', yscalelog=False)
    axs[3,2].set_ylim([0, 1e16])

    for i in range(axs.shape[0]):
        for j in range(axs.shape[1]):
            axs[i,j].set_xlim([t[0], t[-1]])
            axs[i,j].grid(True)


    axs[3, 1].set_ylim([1e-9, 1e-1])


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

def makeplots(so1, so2):
    fig1, axs1 = plt.subplots(4, 3, figsize=(12, 10))

    drawplot1(axs1, so1)
    drawplot1(axs1, so2, toffset=so1.grid.t[-1]+0.7, showlabel=False)

    axs1[0,0].legend(frameon=False)
    axs1[0,1].legend(frameon=False)
    axs1[1,0].legend(frameon=False)
    axs1[2,0].legend(frameon=False)
    axs1[2,1].legend(frameon=False)
    axs1[2,2].legend(frameon=False)
    axs1[3,2].legend(frameon=False)

    fig1.tight_layout()
    plt.show()

def savePlots(so1, so2, directory, filename):
    fig1, axs1 = plt.subplots(4, 3, figsize=(12, 10))

    drawplot1(axs1, so1)
    drawplot1(axs1, so2, toffset=so1.grid.t[-1] + 0.7, showlabel=False)

    axs1[0, 0].legend(frameon=False)
    axs1[0, 1].legend(frameon=False)
    axs1[1, 0].legend(frameon=False)
    axs1[2, 0].legend(frameon=False)
    axs1[2, 1].legend(frameon=False)
    axs1[2, 2].legend(frameon=False)
    axs1[3, 2].legend(frameon=False)

    fig1.tight_layout()
    fig1.savefig(directory + '/' + filename + '.pdf')
    plt.close()

def parameterSweepSingle(n0_list=np.array([]), Vloop_list=np.array([]), Btor_list=np.array([])):
    directory = '../../../Figures/RunawayParameterSweep/Y_Fe^D=1e-5'
    for tritium, nt in zip([False], [3e4]):
        addT = '_T' if tritium else ''
        for n0 in n0_list:
            ss1 = generate(n0=n0, EfieldDyon=False, tritium=tritium)
            ss1.save(f'Sweep/settings1_n0_{np.round(n0,-14)}' + addT + '.h5')
            #ss1.solver.setVerbose(True)
            so1 = runiface(ss1, f'Sweep/output1_n0_{np.round(n0,-14)}' + addT + '.h5', quiet=False)

            ss2 = STREAMSettings(ss1)
            ss2.fromOutput(f'Sweep/output1_n0_{np.round(n0,-14)}' + addT + '.h5')
            ss2.timestep.setTmax(1.3*6 - ss1.timestep.tmax)
            ss2.timestep.setNumberOfSaveSteps(0)
            ss2.timestep.setNt(nt)
            ss2.save(f'Sweep/settings2_n0_{np.round(n0,-14)}' + addT + '.h5')
            so2 = runiface(ss2, f'Sweep/output2_n0_{np.round(n0,-14)}' + addT + '.h5', quiet=False)

            filename = f'n0_{np.round(n0,-14)}'
            savePlots(so1, so2, directory, filename)

        for Vloop in Vloop_list:
            ss1 = generate(Vloop=Vloop, EfieldDyon=False, tritium=tritium)
            ss1.save(f'Sweep/settings1_Vloop_{np.round(Vloop,3)}' + addT + '.h5')
            so1 = runiface(ss1, f'Sweep/output1_Vloop_{np.round(Vloop,3)}' + addT + '.h5', quiet=False)

            ss2 = STREAMSettings(ss1)
            ss2.fromOutput(f'Sweep/output1_Vloop_{np.round(Vloop,3)}' + addT + '.h5')
            ss2.timestep.setTmax(1.3*6 - ss1.timestep.tmax)
            ss2.timestep.setNumberOfSaveSteps(0)
            ss2.timestep.setNt(nt*3)
            ss2.save(f'Sweep/settings2_Vloop_{np.round(Vloop,3)}' + addT + '.h5')
            so2 = runiface(ss2, f'Sweep/output2_Vloop_{np.round(Vloop,3)}' + addT + '.h5', quiet=False)

            filename = f'Vloop_{np.round(Vloop,3)}'
            savePlots(so1, so2, directory, filename)
            
        for Btor in Btor_list:
            ss1 = generate(Btor=Btor, EfieldDyon=False, tritium=tritium)
            ss1.save(f'Sweep/settings1_Btor{np.round(Btor,3)}' + addT + '.h5')
            so1 = runiface(ss1, f'Sweep/output1_Btor{np.round(Btor,3)}' + addT + '.h5', quiet=False)

            ss2 = STREAMSettings(ss1)
            ss2.fromOutput(f'Sweep/output1_Btor{np.round(Btor,3)}' + addT + '.h5')
            ss2.timestep.setTmax(1.3*6 - ss1.timestep.tmax)
            ss2.timestep.setNumberOfSaveSteps(0)
            ss2.timestep.setNt(30000)
            ss2.save(f'Sweep/settings2_Btor{np.round(Btor,3)}' + addT + '.h5')
            so2 = runiface(ss2, f'Sweep/output2_Btor{np.round(Btor,3)}' + addT + '.h5', quiet=False)

            filename = f'Btor{np.round(Btor,3)}'
            savePlots(so1, so2, directory, filename)


def parameterSweepDouble(n0_list=np.array([]), Vloop_list=np.array([]), Btor_list=np.array([])):
    directory = '../../../Figures/RunawayParameterSweep'
    for tritium, nt in zip([False], [1e5]):
        addT = '_T' if tritium else ''
        for Vloop in Vloop_list:
            for n0 in n0_list:
                ss1 = generate(n0=n0, Vloop=Vloop, EfieldDyon=False, tritium=tritium)
                ss1.save(f'Sweep/settings1_n0_{np.round(n0, -14)}_Vloop_{np.round(Vloop, 3)}' + addT + '.h5')
                so1 = runiface(ss1, f'Sweep/output1_n0_{np.round(n0, -14)}_Vloop_{np.round(Vloop, 3)}' + addT + '.h5', quiet=False)

                ss2 = STREAMSettings(ss1)
                ss2.fromOutput(f'Sweep/output1_n0_{np.round(n0, -14)}_Vloop_{np.round(Vloop, 3)}' + addT + '.h5')
                ss2.timestep.setTmax(1.3 * 6 - ss1.timestep.tmax)
                ss2.timestep.setNumberOfSaveSteps(0)
                ss2.timestep.setNt(nt)
                ss2.save(f'Sweep/settings2_n0_{np.round(n0, -14)}_Vloop_{np.round(Vloop, 3)}' + addT + '.h5')
                so2 = runiface(ss2, f'Sweep/output2_n0_{np.round(n0, -14)}_Vloop_{np.round(Vloop, 3)}' + addT + '.h5', quiet=False)

                filename = f'n0_{np.round(n0, -14)}'
                savePlots(so1, so2, directory, filename)

def oneRun(argv):
    FONTSIZE = 16
    plt.rcParams.update({'font.size': FONTSIZE})

    parser = argparse.ArgumentParser()

    parser.add_argument('-e', '--extension', help="Append the specified string to the end of file names", dest="extension", action='store', default='')
    parser.add_argument('-n', '--no-plot', help="Only run simulations, don't generate plots", dest="plot", action='store_false', default=True)
    parser.add_argument('-s', '--skip', help="Skip the simulation and load output from already existing files", action='store_true', dest="skip", default=False)

    settings = parser.parse_args()

    ext = '' if not settings.extension else '_' + settings.extension

    if not settings.skip:
        ss1 = generate(n0=1.5e17, EfieldDyon=False, tritium=False)
        ss1.save(f'settings1{ext}.h5')
        so1 = runiface(ss1, f'output1{ext}.h5', quiet=False)
        print('ye')
        ss2 = STREAMSettings(ss1)
        ss2.fromOutput(f'output1{ext}.h5')
        ss2.timestep.setTmax(1.3*6 - ss1.timestep.tmax)
        ss2.timestep.setNumberOfSaveSteps(0)
        ss2.timestep.setNt(30000)
        ss2.save(f'settings2{ext}.h5')
        so2 = runiface(ss2, f'output2{ext}.h5', quiet=False)
    else:
        so1 = STREAMOutput(f'output1{ext}.h5')
        so2 = STREAMOutput(f'output2{ext}.h5')

    if settings.plot:
        makeplots(so1, so2)

    return 0

def main(argv):
    oneRun(argv=argv)
    '''
    n0_list = np.array([1.5, 1.75]) * 1e17
    #Vloop_list = np.array([2, 3])
    #Btor_list = np.array([2.5, 3])
    parameterSweepSingle(n0_list=n0_list)#, Vloop_list=Vloop_list)
    #'''

    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


