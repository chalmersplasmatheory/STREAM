#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import time
import scipy.constants
import sys

sys.path.append('../../../py')

from PYDYON import Simulation
from STREAM import STREAMSettings, STREAMOutput

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

OUTFILE = '2'

ss = STREAMSettings(f'settings{OUTFILE}.h5')
so = STREAMOutput(f'output{OUTFILE}.h5')


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
    'Lp': scipy.constants.mu_0*R0*(np.log(8*R0/a[0]) + 0.5 - 2),
    'LMK2': 9.1e-6,
    'M': 2.49e-6,
    'RMK2': 7.5e-4,
    'simple': False
}

print('PYDYON: \n'+str(params)+'\n\n')
sparam = loadSTREAMSettings(ss)
print('STREAM: \n'+str(sparam)+'\n\n')

sim = Simulation(**sparam)

sim.addIon('D', 1)
sim.addIon('C', 6)
sim.addIon('O', 8)
sim.initialize(Te=1, Ti=0.03, Ip=2.4e3, IMK2=0, niD=nD, niC=nC, niO=nO)

tic = time.time()
solution = sim.solve(tMax=0.3)
#solution = sim.solve(tMax=1e-6)

print('Obtained solution in {:.3f} s'.format(time.time()-tic))

#print(list(solution.keys()))
#solution.plotIonDensity('D')
#solution.plotIonDensity('C')
#solution.plotIonDensity('O')
plt.show()
solution.plot()
#solution.plotJET()
#solution.plotKimThesis45()
solution.plotKim2020()
#'''
