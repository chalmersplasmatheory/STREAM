#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import time
import scipy.constants
import sys

sys.path.append('../py')

from PYDYON import ConfinementTime, IonHandler, PlasmaVolume, UnknownQuantityHandler, Simulation, SimulationResult
from STREAM import STREAMSettings, STREAMOutput
from PYDYON.Equations import *

SETTINGS = {
    'a': lambda ss : ss.radialgrid.a,
    'ta': lambda ss : ss.radialgrid.ta,
    'R': lambda ss : ss.radialgrid.R0,
    'V_vessel': lambda ss : ss.radialgrid.vessel_volume,
    #'kappa': lambda ss : ss.radialgrid.kappa,
    'Bphi': lambda ss : ss.radialgrid.B0[0],
    'Bv': lambda ss : 1e-3,
    # The following currently only work with the 'TYPE_CIRCUIT' model in STREAM
    'l_MK2': lambda ss : ss.radialgrid.b,
    'Vloop': lambda ss : ss.eqsys.E_field.circuit_Vloop if ss.eqsys.E_field.circuit_Vloop.size==1 else scipy.interpolate.interp1d(ss.eqsys.E_field.circuit_Vloop_t, ss.eqsys.E_field.circuit_Vloop, bounds_error = False, fill_value = 'extrapolate'),
    'Lp': lambda ss : ss.eqsys.E_field.circuit_Lp,
    'LMK2': lambda ss : ss.eqsys.E_field.circuit_Lwall,
    'M': lambda ss : ss.eqsys.E_field.circuit_M,
    'RMK2': lambda ss : ss.eqsys.E_field.circuit_Rwall if ss.eqsys.E_field.circuit_Rwall<1 else np.inf,
    # "simple" neutral deuterium influx
    #S'simple': lambda ss : ss.radialgrid.c2==1
}
def loadSTREAMSettings(ss):
    """
    Convert a STREAMSettings object into a dict with settings
    for PYDYON.
    """
    # Load settings from STREAMSettings object
    settings = {}
    for s in SETTINGS.keys():
        settings[s] = SETTINGS[s](ss)

    return settings

ss = STREAMSettings('settings_final.h5')
so = STREAMOutput('output_final.h5')

settings = loadSTREAMSettings(ss)

sim = Simulation(**settings)

prefill = 2.7e-3    # Pa
gamma_i = 2e-3      # Ionization fraction
nD0 = 4.8e20 * prefill
nD = nD0 * np.array([1 - gamma_i, gamma_i])

nC = np.zeros((7,))
nO = np.zeros((9,))
nO[0] = 1e-3 * nD0

# Construct ion object
for ion in so.eqsys.n_i.ions:
    sim.addIon(ion.name,ion.Z)

sim.initialize(Te=1, Ti=0.03, Ip=2.4e3, IMK2=0, niD=nD, niC=nC, niO=nO)

solution = sim.solve(tMax=0.3)
uqh = solution.simulation.unknowns

RTOL = 1e-1

'''
# lambda_D
name = 'lambda_D'
stream = so.eqsys.lambda_i['D'][1:,0]
pv     = PlasmaVolume(a=settings['a'], R=settings['R'], V_vessel=settings['V_vessel'], ions=ions, ta=settings['ta'])
pydyon = pv.getLambda('D', uqh['ne'], uqh['Te'], uqh['Ti'])
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(so.grid.t[1:], stream, 'k', label='STREAM')
ax.plot(solution.getT(), pydyon, 'r--', label='PYDYON')
# ax.plot(t, stream.flatten()/pydyon.flatten(), 'k')
ax.set_xlabel('Time (s)')
ax.set_title(name)
ax.set_xlim([0, 0.3])
ax.legend()
fig.tight_layout()
plt.show()

# lambda_C
name = 'lambda_C'
stream = so.eqsys.lambda_i['C'][1:,0]
pv     = PlasmaVolume(a=settings['a'], R=settings['R'], V_vessel=settings['V_vessel'], ions=ions, ta=settings['ta'])
pydyon = pv.getLambda('C', uqh['ne'], uqh['Te'], uqh['Ti'])
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(so.grid.t[1:], stream, 'k', label='STREAM')
ax.plot(solution.getT(), pydyon, 'r--', label='PYDYON')
# ax.plot(t, stream.flatten()/pydyon.flatten(), 'k')
ax.set_xlabel('Time (s)')
ax.set_title(name)
ax.set_xlim([0, 0.3])
ax.legend()
fig.tight_layout()
plt.show()

# lambda_O
name = 'lambda_O'
stream = so.eqsys.lambda_i['O'][1:,0]
pv     = PlasmaVolume(a=settings['a'], R=settings['R'], V_vessel=settings['V_vessel'], ions=ions, ta=settings['ta'])
pydyon = pv.getLambda('O', uqh['ne'], uqh['Te'], uqh['Ti'])
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(so.grid.t[1:], stream, 'k', label='STREAM')
ax.plot(solution.getT(), pydyon, 'r--', label='PYDYON')
# ax.plot(t, stream.flatten()/pydyon.flatten(), 'k')
ax.set_xlabel('Time (s)')
ax.set_title(name)
ax.set_xlim([0, 0.3])
ax.legend()
fig.tight_layout()
plt.show()
'''
i = 1
# Radiated power
name = 'Radiated power'
stream = so.other.fluid.Tcold_radiation[:,0]
pterm  = RadiatedPowerTerm(sim.unknowns, sim.ions)
pydyon = solution.evaluateTerm(pterm)
fig = plt.figure(i)
ax = fig.add_subplot(111)
ax.plot(so.grid.t[1:], stream, 'k', label='STREAM')
ax.plot(solution.getT(), pydyon, 'r--', label='PYDYON')
# ax.plot(t, stream.flatten()/pydyon.flatten(), 'k')
ax.set_xlabel('Time (s)')
ax.set_title(name)
ax.set_xlim([0, 0.3])
ax.legend()
fig.tight_layout()
i=i+1

# Ohmic Power
name = 'Ohmic power'
stream = -so.other.fluid.Tcold_ohmic[:,0]
pterm  = OhmicPowerTerm(sim.unknowns, sim.ions)
pydyon = solution.evaluateTerm(pterm)
fig = plt.figure(i)
ax = fig.add_subplot(111)
ax.plot(so.grid.t[1:], stream, 'k', label='STREAM')
ax.plot(solution.getT(), pydyon, 'r--', label='PYDYON')
# ax.plot(t, stream.flatten()/pydyon.flatten(), 'k')
ax.set_xlabel('Time (s)')
ax.set_title(name)
ax.set_xlim([0, 0.3])
ax.legend()
fig.tight_layout()
i=i+1

# e-i equilibration
name = 'e-i equilibration'
stream = so.other.fluid.Tcold_ion_coll[:,0]
pterm  = EquilibrationPowerTerm(sim.unknowns, sim.ions)
pydyon = solution.evaluateTerm(pterm)
fig = plt.figure(i)
ax = fig.add_subplot(111)
ax.plot(so.grid.t[1:], stream, 'k', label='STREAM')
ax.plot(solution.getT(), pydyon, 'r--', label='PYDYON')
# ax.plot(t, stream.flatten()/pydyon.flatten(), 'k')
ax.set_xlabel('Time (s)')
ax.set_title(name)
ax.set_xlim([0, 0.3])
ax.legend()
fig.tight_layout()
i=i+1

# e heat convection
name = 'e heat convection'
stream = so.other.stream.Tcold_transport[:,0]
pterm  = ElectronConvectivePowerTerm(sim.unknowns, sim.ions)
pydyon = solution.evaluateTerm(pterm)
fig = plt.figure(i)
ax = fig.add_subplot(111)
ax.plot(so.grid.t[1:], stream, 'k', label='STREAM')
ax.plot(solution.getT(), pydyon, 'r--', label='PYDYON')
# ax.plot(t, stream.flatten()/pydyon.flatten(), 'k')
ax.set_xlabel('Time (s)')
ax.set_title(name)
ax.set_xlim([0, 0.3])
ax.legend()
fig.tight_layout()
i=i+1

# i heat convection
name = 'i heat convection'
stream = -np.sum(so.other.stream.Wi_iontransport[:,:,0], axis=1)
pterm  = IonConvectivePowerTerm(sim.unknowns, sim.ions)
pydyon = solution.evaluateTerm(pterm)
fig = plt.figure(i)
ax = fig.add_subplot(111)
ax.plot(so.grid.t[1:], stream, 'k', label='STREAM')
ax.plot(solution.getT(), pydyon, 'r--', label='PYDYON')
# ax.plot(t, stream.flatten()/pydyon.flatten(), 'k')
ax.set_xlabel('Time (s)')
ax.set_title(name)
ax.set_xlim([0, 0.3])
ax.legend()
fig.tight_layout()
i=i+1

# Charge-exchange heat loss
name = 'Charge-exchange heat loss'
stream = -so.other.stream.Wi_chargeexchange[:,0]
pterm  = ChargeExchangePowerTerm(sim.unknowns, sim.ions)
pydyon = solution.evaluateTerm(pterm)
fig = plt.figure(i)
ax = fig.add_subplot(111)
ax.plot(so.grid.t[1:], stream, 'k', label='STREAM')
ax.plot(solution.getT(), pydyon, 'r--', label='PYDYON')
# ax.plot(t, stream.flatten()/pydyon.flatten(), 'k')
ax.set_xlabel('Time (s)')
ax.set_title(name)
ax.set_xlim([0, 0.3])
ax.legend()
fig.tight_layout()
i=i+1
'''
# i particle transport
name = 'i particle transport'
stream = -so.other.stream.ni_iontransport[:,1,0]
pterm  = IonTransport(sim.unknowns, sim.ions)
pydyon = solution.evaluateTerm(pterm)
fig = plt.figure(i)
ax = fig.add_subplot(111)
ax.plot(so.grid.t[1:], stream, 'k', label='STREAM')
ax.plot(solution.getT(), pydyon, 'r--', label='PYDYON')
# ax.plot(t, stream.flatten()/pydyon.flatten(), 'k')
ax.set_xlabel('Time (s)')
ax.set_title(name)
ax.set_xlim([0, 0.3])
ax.legend()
fig.tight_layout()
i=i+1
'''
# Confinement time
name = 'Confinement time'
stream = so.other.stream.tau_D[:,0]
pterm  = ConfinementTime(sim.unknowns, sim.ions)
pydyon = solution.evaluateTerm(pterm)
fig = plt.figure(i)
ax = fig.add_subplot(111)
ax.plot(so.grid.t[1:], stream, 'k', label='STREAM')
ax.plot(solution.getT(), pydyon, 'r--', label='PYDYON')
# ax.plot(t, stream.flatten()/pydyon.flatten(), 'k')
ax.set_xlabel('Time (s)')
ax.set_title(name)
ax.set_xlim([0, 0.3])
ax.legend()
fig.tight_layout()
i=i+1




plt.show()