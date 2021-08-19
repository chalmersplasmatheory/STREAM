#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import time
import sys

sys.path.append('../../py')

from PYDYON import Simulation


prefill = 5e-5
gamma_i = 2e-3      # Ionization fraction
nD0 = 3.22e22 * prefill
nD = nD0 * np.array([1-gamma_i, gamma_i])

sim = Simulation(RMK2=1e5)

sim.addIon('D', 1)
sim.initialize(Te=1, Ti=0.03, Ip=450, IMK2=0, niD=nD)

tic = time.time()
#solution = sim.solve(tMax=0.03)
solution = sim.solve(tMax=1e-6)

print('Obtained solution in {:.3f} s'.format(time.time()-tic))

#print(list(solution.keys()))
solution.plot()

