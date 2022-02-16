# Generate baseline STREAM settings

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import mu_0
import sys

#sys.path.append('/mnt/HDD/software/STREAM/py')
#sys.path.append('/run/media/mathias/BlackTools/Fusion/software/STREAM/py')
sys.path.append('../../py')

import PlasmaParameters as Formulas
from STREAM import STREAMSettings, runiface
import DREAM.Settings.Atomics as Atomics
import STREAM.Settings.Equations.ElectricField as ElectricField
import STREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.Equations.ColdElectronTemperature as Tcold


def generate_baseline(prefill=1e-6, tritium=True, Vloop0=12, runaways=True, Ip0=40e3, gamma=2e-3):
    """
    Generate an ITER baseline scenario.

    :param prefill:       Prefill gas pressure [Torr].
    :param bool tritium:  Whether or not to include tritium in the simulation.
    :param Vloop:         Externally applied loop voltage.
    :param bool runaways: Whether or not to include runaway generation in the simulation.
    :param Ip0:           Initial plasma current.
    :param gamma:         Initial plasma ionization fraction.
    """
    ss = STREAMSettings()

    n0 = 3.22e22 * prefill

    if tritium:
        nD = np.array([[n0*0.5], [n0 * gamma / (1-gamma)*0.5]])
        nT = np.array([[n0*0.5], [n0 * gamma / (1-gamma)*0.5]])
    else:
        nD = np.array([[n0], [n0 * gamma / (1-gamma)]])
        nT = 0

    a = 1.6                 # Plasma minor radius [m]
    R0 = 5.65               # Plasma major radius [m]
    #b = 8*R0/np.exp(7/4)    # TODO
    b = 2
    V_vessel = 1700         # Vacuum vessel volume [m^3]
    Btor = 2.65             # Toroidal magnetic field strength [T]

    Te0 = 5     # initial electron temperature [eV]
    Ti0 = 1     # initial ion temperature [eV]

    j0 = Ip0 / (a**2 * np.pi)
    E0 = j0 / Formulas.evaluateSpitzerConductivity(n=nD[1], T=Te0, Z=1)

    ss.atomic.adas_interpolation = Atomics.ADAS_INTERP_BILINEAR

    # Electric field dynamics
    #Lp = float(mu_0 * R0 * (np.log(8*R0/a) + 0.25 - 2))
    #print(f'Lp = {Lp},  DREAM Lp = {mu_0*R0}')
    Lp = mu_0*R0
    Lwall = 9.1e-6
    Rwall = 10e6
    M = 2.49e-6

    ss.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
    ss.eqsys.E_field.setInitialProfile(E0)
    ss.eqsys.E_field.setInductances(Lp=Lp, Lwall=Lwall, M=M, Rwall=Rwall)
    ss.eqsys.E_field.setCircuitVloop(Vloop0)

    # Electron temperature
    ss.eqsys.T_cold.setType(Tcold.TYPE_SELFCONSISTENT)
    ss.eqsys.T_cold.setInitialProfile(Te0)
    ss.eqsys.T_cold.setRecombinationRadiation(Tcold.RECOMBINATION_RADIATION_INCLUDED)

    # Ions
    ss.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=nD, r=np.array([0]), T=Ti0)
    if tritium:
        ss.eqsys.n_i.addIon(name='T', Z=1, iontype=Ions.IONS_DYNAMIC, n=nT, r=np.array([0]), T=Ti0)

    # Recycling coefficients

    # Enable runaway?
    if runaways:
        ss.eqsys.n_re.setAvalanche(Runaways.AVALANCHE_MODE_FLUID_HESSLOW)
        ss.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_NEURAL_NETWORK)

        if tritium:
            ss.eqsys.n_re.setTritium(True)
            ss.eqsys.n_re.setCompton(Runaways.COMPTON_MODE_FLUID, Runaways.ITER_PHOTON_FLUX_DENSITY)

    # Radial grid
    ss.radialgrid.setB0(Btor)
    ss.radialgrid.setMinorRadius(a)
    ss.radialgrid.setMajorRadius(R0)
    ss.radialgrid.setWallRadius(b)
    ss.radialgrid.setVesselVolume(V_vessel)

    ss.radialgrid.setRecyclingCoefficient1(1)
    ss.radialgrid.setRecyclingCoefficient2(0)
    ss.radialgrid.setRecyclingCoefficient3(1)

    # Disable kinetic grids
    ss.hottailgrid.setEnabled(False)
    ss.runawaygrid.setEnabled(False)

    ss.solver.setType(Solver.NONLINEAR)
    ss.timestep.setTmax(1e-4)
    ss.timestep.setNt(10000)
    ss.timestep.setNumberOfSaveSteps(10000)

    ss.other.include('fluid', 'stream', 'scalar')

    return ss


