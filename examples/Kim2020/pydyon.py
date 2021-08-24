#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import time
import scipy.constants
import sys

sys.path.append('../../py')

from PYDYON import Simulation


prefill = 2 * 1.2e-3 / 133.32
print('Prefill pressure: {:.2e} Torr'.format(prefill))
gamma_i = 2e-3      # Ionization fraction
#nD0 = 3.22e22 * prefill
nD0 = 2.78e22 * prefill
nD = nD0 * np.array([1-gamma_i, gamma_i])

params = {
    'a': 1.6,
    'R': 5.65,
    'V_vessel': 1000,
    'Bphi': 2.65,
    'Bv': 2e-3,
    'l_MK2': 2,
    'Vloop': 12,
    'Lp': scipy.constants.mu_0*5.65*(np.log(8*5.65/1.6) + 0.5 - 2),
    'RMK2': np.inf
}

sim = Simulation(**params)

sim.addIon('D', 1)
sim.initialize(Te=1, Ti=0.03, Ip=2.4e3, IMK2=0, niD=nD)

tic = time.time()
solution = sim.solve(tMax=0.3)

print('Obtained solution in {:.3f} s'.format(time.time()-tic))

#print(list(solution.keys()))
solution.plot()
#solution.plotKimThesis45()
solution.plotKim2020()

