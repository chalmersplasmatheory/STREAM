#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import time
import scipy.constants
from scipy.interpolate import interp1d
import sys

sys.path.append('../py')

from PYDYON import ConfinementTime, IonHandler, PlasmaVolume, UnknownQuantityHandler, Simulation, SimulationResult, STREAM
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
    #'Vloop': lambda ss : ss.eqsys.E_field.circuit_Vloop if ss.eqsys.E_field.circuit_Vloop.size==1 else scipy.interpolate.interp1d(ss.eqsys.E_field.circuit_Vloop_t, ss.eqsys.E_field.circuit_Vloop, bounds_error = False, fill_value = 'extrapolate'),
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

#ss = STREAMSettings('settings_final.h5')
so = STREAMOutput('output_final.h5')

#settings = loadSTREAMSettings(ss)

prefill = 2.7e-3    # Pa
gamma_i = 2e-3      # Ionization fraction
nD0 = 4.8e20 * prefill
nD = nD0 * np.array([1-gamma_i, gamma_i])

nC = np.zeros((7,))
nO = np.zeros((9,))
nO[0] = 1e-3 * nD0

R0 = 2.96

t_a   = [0  ,  0.017,  0.05,  0.085,  0.14,  0.19,  0.25,  0.3]
V_p   = np.array([100, 80    , 56   , 48    , 52   , 51.75, 54.25, 56  ])
tet   = np.linspace(0, 0.3, 100)
a_fun   = scipy.interpolate.interp1d(t_a, np.sqrt(V_p/(2*np.pi**2*R0)))
a=a_fun(tet)

t_Vloop = [0 , 0.02 , 0.0325, 0.0475, 0.08, 0.1 , 0.125, 0.13, 0.15, 0.20, 0.22, 0.23, 0.25, 0.3 , 0.335, 0.35, 0.37, 0.4 , 0.45, 0.5 ]
d_Vloop = [11, 21.25, 26    , 26.25 , 24  , 16.5, 8.25 , 7.9 , 7.75, 7.5 , 7.25, 6.5 , 6.5 , 6.75, 6.75 , 6   , 4.75, 4.25, 4.5 , 3.60]
Vloop   = scipy.interpolate.interp1d(t_Vloop, d_Vloop)

params = {
    'a': a,
    'ta': tet,
    'R': R0,
    'V_vessel': 100,
    'Bphi': 2.4,
    'Bv': 1e-3,
    'l_MK2': 1,
    'Vloop': Vloop,
    'Lp': 5.19e-6,
    'LMK2': 9.1e-6,
    'M': 2.49e-6,
    'RMK2': 7.5e-4
}

#sim = Simulation(**settings)
sim = Simulation(**params)

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
ions = solution.simulation.ions

RTOL = 1e-1
i = 1

solution.plotJetPydyonVsStream(so)

'''
# Lambda_O
name = 'Lambda_O'
stream = so.eqsys.lambda_i['O'][1:,0]
pydyon = uqh.plasmavolume.getLambdaVec('O', solution.x['ne'][:], solution.x['Te'][:], solution.x['Ti'][:])
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
'''
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
plt.show() #i=i+1

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
plt.show() #i=i+1
#'''
'''
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
plt.show() #i=i+1

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
plt.show() #i=i+1

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
plt.show() #i=i+1

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
plt.show() #i=i+1
#'''
'''
# D particle transport
name = 'i particle transport'
stream = -so.other.stream.ni_iontransport[:,1,0]
pterm  = IonTransport(sim.unknowns, sim.ions)
pydyon = solution.evaluateTermIonState(pterm, 'D', 1)
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
plt.show() #i=i+1
#'''
'''
# O particle transport
name = 'O particle transport'
stream = 0
pterm  = IonTransport(sim.unknowns, sim.ions)
pydyon = 0
for i in range(1,9):
    stream = stream-so.other.stream.ni_iontransport[:,i+9,0]
    pydyon = pydyon+solution.evaluateTermIonState(pterm, 'O', i)
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
#'''
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
plt.show() #i=i+1
#'''
'''
# D0 density
name = 'D0 density'
stream = so.eqsys.n_i['D'][0][1:]
pydyon = solution.x['niD'][0,:]
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
plt.show() #i=i+1

# D1 density
name = 'D1 density'
stream = so.eqsys.n_i['D'][1][1:]
pydyon = solution.x['niD'][1,:]
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
plt.show() #i=i+1
#'''
'''
# C0 density
name = 'C0 density'
stream = so.eqsys.n_i['C'][0][1:]
pydyon = solution.x['niC'][0,:]
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
plt.show() #i=i+1

# C1 density
name = 'C1 density'
stream = so.eqsys.n_i['C'][1][1:]
pydyon = solution.x['niC'][1,:]
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
plt.show() #i=i+1

# C2 density
name = 'C2 density'
stream = so.eqsys.n_i['C'][2][1:]
pydyon = solution.x['niC'][2,:]
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
plt.show() #i=i+1

# C3 density
name = 'C3 density'
stream = so.eqsys.n_i['C'][3][1:]
pydyon = solution.x['niC'][3,:]
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
plt.show() #i=i+1

# C4 density
name = 'C4 density'
stream = so.eqsys.n_i['C'][4][1:]
pydyon = solution.x['niC'][4,:]
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
plt.show() #i=i+1

# C5 density
name = 'C5 density'
stream = so.eqsys.n_i['C'][5][1:]
pydyon = solution.x['niC'][5,:]
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
plt.show() #i=i+1

# C6 density
name = 'C6 density'
stream = so.eqsys.n_i['C'][6][1:]
pydyon = solution.x['niC'][6,:]
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
plt.show() #i=i+1

# O0 density
name = 'O0 density'
stream = so.eqsys.n_i['O'][0][1:]
pydyon = solution.x['niO'][0,:]
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
plt.show() #i=i+1

# O1 density
name = 'O1 density'
stream = so.eqsys.n_i['O'][1][1:]
pydyon = solution.x['niO'][1,:]
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
plt.show() #i=i+1

# O2 density
name = 'O2 density'
stream = so.eqsys.n_i['O'][2][1:]
pydyon = solution.x['niO'][2,:]
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
plt.show() #i=i+1

# O3 density
name = 'O3 density'
stream = so.eqsys.n_i['O'][3][1:]
pydyon = solution.x['niO'][3,:]
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
plt.show() #i=i+1

# O4 density
name = 'O4 density'
stream = so.eqsys.n_i['O'][4][1:]
pydyon = solution.x['niO'][4,:]
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
plt.show() #i=i+1

# O5 density
name = 'O5 density'
stream = so.eqsys.n_i['O'][5][1:]
pydyon = solution.x['niO'][5,:]
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
plt.show() #i=i+1

# O6 density
name = 'O6 density'
stream = so.eqsys.n_i['O'][6][1:]
pydyon = solution.x['niO'][6,:]
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
plt.show() #i=i+1

# O7 density
name = 'O7 density'
stream = so.eqsys.n_i['O'][7][1:]
pydyon = solution.x['niO'][7,:]
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
plt.show() #i=i+1

# O8 density
name = 'O8 density'
stream = so.eqsys.n_i['O'][8][1:]
pydyon = solution.x['niO'][8,:]
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
plt.show() #i=i+1
#'''
'''
# D0 negative ionization
name = 'D0 negative ionization'
stream = so.other.stream.ionrateequation_negIonization[:,0,0]
pterm  = DeuteriumAtomBalance(sim.unknowns, sim.ions)
pydyon = solution.evaluateTermFull(pterm, 2, True)
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
plt.show() #i=i+1

# D0 positive recombination
name = 'D0 positive recombination'
stream = so.other.stream.ionrateequation_posRecombination[:,0,0]
pterm  = DeuteriumAtomBalance(sim.unknowns, sim.ions)
pydyon = solution.evaluateTermFull(pterm, 3, True)
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
plt.show() #i=i+1
#'''
'''
# D0 negative C-X
name = 'D0 negative C-X'
stream = so.other.stream.ionrateequation_negChargeExchange[:,0,0]
pterm  = DeuteriumAtomBalance(sim.unknowns, sim.ions)
pydyon = solution.evaluateTermFull(pterm, 6, True)
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
plt.show() #i=i+1
#'''
'''
# D0 influx
name = 'D0 influx'
stream = so.other.stream.neutralinflux['D'][:] / so.other.stream.V_n_tot['D'][:]
pterm  = DeuteriumInflux(sim.unknowns, sim.ions)
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
plt.show() #i=i+1
#'''
'''
# O0 influx
name = 'O0 influx'
stream = so.other.stream.neutralinflux['O'][:] / so.other.stream.V_n_tot['O'][:]
pterm  = IonInflux(sim.unknowns, sim.ions, )
pydyon = solution.evaluateTermIon(pterm, 'O')
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
plt.show() #i=i+1
#'''
'''
# D1 positive ionization
name = 'D1 positive ionization'
stream = so.other.stream.ionrateequation_posIonization[:,1,0]
pterm  = DeuteriumIonBalance(sim.unknowns, sim.ions)
pydyon = solution.evaluateTermFull(pterm, 1, True)
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
plt.show() #i=i+1


# D1 negative recombination
name = 'D1 negative recombination'
stream = so.other.stream.ionrateequation_negRecombination[:,1,0]
pterm  = DeuteriumIonBalance(sim.unknowns, sim.ions)
pydyon = solution.evaluateTermFull(pterm, 4, True)
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
plt.show() #i=i+1
#'''
'''
# D1 positive C-X
name = 'D1 positive C-X'
stream = so.other.stream.ionrateequation_posChargeExchange[:,1,0]
pterm  = DeuteriumIonBalance(sim.unknowns, sim.ions)
pydyon = solution.evaluateTermFull(pterm, 5, True)
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
plt.show() #i=i+1
#'''

#plt.show()

'''
# D0 negative recombination
name = 'D0 negative recombination'
stream = so.other.stream.ionrateequation_negRecombination[:,0,0]
pterm  = DeuteriumAtomBalance(sim.unknowns, sim.ions)
pydyon = solution.evaluateTermFull(pterm, 4, True)
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

# D0 positive C-X
name = 'D0 positive C-X'
stream = so.other.stream.ionrateequation_posChargeExchange[:,0,0]
pterm  = DeuteriumAtomBalance(sim.unknowns, sim.ions)
pydyon = solution.evaluateTermFull(pterm, 5, True)
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

# D1 negative ionization
name = 'D1 negative ionization'
stream = so.other.stream.ionrateequation_negIonization[:,1,0]
pterm  = DeuteriumIonBalance(sim.unknowns, sim.ions)
pydyon = solution.evaluateTermFull(pterm, 2, True)
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

# D1 positive recombination
name = 'D1 positive recombination'
stream = so.other.stream.ionrateequation_posRecombination[:,1,0]
pterm  = DeuteriumIonBalance(sim.unknowns, sim.ions)
pydyon = solution.evaluateTermFull(pterm, 3, True)
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

# D1 negative C-X
name = 'D1 negative C-X'
stream = so.other.stream.ionrateequation_negChargeExchange[:,1,0]
pterm  = DeuteriumIonBalance(sim.unknowns, sim.ions)
pydyon = solution.evaluateTermFull(pterm, 6, True)
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

# D0 positive ionization
name = 'D0 positive ionization'
stream = so.other.stream.ionrateequation_posIonization[:,0,0]
pterm  = DeuteriumAtomBalance(sim.unknowns, sim.ions)
pydyon = solution.evaluateTermFull(pterm, 1, True)
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