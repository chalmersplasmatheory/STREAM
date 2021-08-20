#!/usr/bin/python3

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
opt = PYDYON.Equations.OhmicPowerTerm(uqh, ions)
ept = PYDYON.Equations.EquilibrationPowerTerm(uqh, ions)
ecpt = PYDYON.Equations.ElectronConvectivePowerTerm(uqh, ions)
tau = PYDYON.ConfinementTime(uqh, ions)
cir = PYDYON.Equations.CircuitEquation(uqh, ions, Vloop=20, RMK2=1e6)

so = STREAMOutput('output12.h5')
#so = STREAMOutput('output11.h5')

# Evaluate Prad
t = so.grid.t[1:]
P = np.zeros((t.size,))
for i in range(t.size):
    x = PYDYON.fromSTREAM(so, uqh, time=i)
    #P[i] = rpt(x)
    #P[i] = opt(x)
    #P[i] = ept(x)
    #P[i] = ecpt(x)
    #P[i] = tau(x)
    P[i] = cir.dIp_dt(t[i], x)

# Radiated power term (suspicious!)
#plt.plot(t, so.other.fluid.Tcold_radiation[:,0], 'k', label='STREAM')
#plt.plot(t, P, 'r--', label='PYDYON')
# Ohmic power term (OK)
#plt.plot(t, -so.other.fluid.Tcold_ohmic[:,0], 'k', label='STREAM')
#plt.plot(t, P, 'r--', label='PYDYON')
# Equilibration power term (OKish)
#plt.plot(t, so.other.fluid.Tcold_ion_coll[:,0], 'k', label='STREAM')
#plt.plot(t, P, 'r--', label='PYDYON')
# Electron convective power term (suspicious!)
#plt.plot(t, so.other.scalar.energyloss_T_cold[:,0], 'k', label='STREAM')
#plt.plot(t, P, 'r--', label='PYDYON')
# Confinement time (OK)
#plt.plot(t, so.other.stream.tau_D[:,0], 'k', label='STREAM')
#plt.plot(t, P, 'r--', label='PYDYON')
# dIp/dt
plt.plot(t, np.diff(so.eqsys.I_p[:,0]) / np.diff(so.grid.t), 'k', label='STREAM')
plt.plot(t, P, 'r--', label='PYDYON')

plt.xlabel(r'Time (s)')
plt.ylabel(r'Power density W/m$^3$')
plt.legend()
plt.tight_layout()

plt.show()

