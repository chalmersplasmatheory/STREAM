#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import time
import sys

sys.path.append('../../py')

from PYDYON import Simulation


#prefill = 3e-5
prefill = 4e-5
#prefill = 7e-5
gamma_i = 2e-3      # Ionization fraction
#nD0 = 3.22e22 * prefill
nD0 = 2.78e22 * prefill
nD = nD0 * np.array([1-gamma_i, gamma_i])

sim = Simulation(RMK2=np.inf, a=0.5)

sim.addIon('D', 1)
sim.initialize(Te=1, Ti=0.03, Ip=450, IMK2=0, niD=nD)

tic = time.time()
solution = sim.solve(tMax=0.1)

print('Obtained solution in {:.3f} s'.format(time.time()-tic))

#print(list(solution.keys()))
solution.plot()
solution.plotKimThesis45()

