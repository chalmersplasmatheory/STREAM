#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../py')

import PYDYON, PYDYON.Equations
from STREAM import STREAMOutput


def setupPYDYON():
    ih  = PYDYON.IonHandler()
    ih.addIon('D', Z=1)

    pv = PYDYON.PlasmaVolume(a=0.5, R=3, V_vessel=100, ions=ih)

    uqh = PYDYON.UnknownQuantityHandler(ih, pv)

    return ih, uqh


ions, uqh = setupPYDYON()

rpt = PYDYON.Equations.RadiatedPowerTerm(uqh, ions)
so = STREAMOutput('output12.h5')

# Evaluate Prad
t = so.grid.t[1:]
P = np.zeros((t.size,))
for i in range(t.size):
    x = PYDYON.fromSTREAM(so, uqh, time=i)
    P[i] = rpt.Prad(x)

plt.plot(t, so.other.fluid.Tcold_radiation[:,0], 'k', label='STREAM')
plt.plot(t, P, 'r--', label='PYDYON')

plt.xlabel(r'Time (s)')
plt.ylabel(r'Power density W/m$^3$')
plt.legend()

plt.show()

