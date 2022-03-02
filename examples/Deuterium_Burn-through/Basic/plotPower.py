#!/usr/bin/env python3
#
# Plot the P_equi collisional power loss term.

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
import sys

sys.path.append('../../py')

from STREAM import STREAMOutput


def combine(fcn, *objs):
    time = []
    data = []

    t0 = 0
    for o in objs:
        time.append(t0 + o.grid.t[1:])
        t0 = time[-1][-1]
        data.append(fcn(o))

    l = sum([t.size for t in time])
    t = np.zeros((l,))
    d = np.zeros((l,))

    offset = 0
    for _t, _d in zip(time, data):
        l = offset+_t.size
        t[offset:l] = _t
        d[offset:l] = _d

        offset = l

    return t, d


def getPequi(*args):
    t, Te = combine(lambda d : d.eqsys.T_cold[1:,0], *args)
    _, Ti = combine(lambda d : d.eqsys.W_i.getTemperature('D')[1:,0], *args)
    _, ne = combine(lambda d : d.eqsys.n_cold[1:,0], *args)
    _, ni = combine(lambda d : d.eqsys.n_i['D'][1][1:,0], *args)
    _, logLambda = combine(lambda d : d.other.fluid.lnLambdaT[:,0], *args)
    M_D = 2     # Deuterium mass in AMU

    P_equi = 7.75e-34 * (Te - Ti) / (Te**(3/2)) * ne*logLambda * ni/M_D

    return t, P_equi


def getTauPar(*args):
    so1 = args[0]

    t, Te = combine(lambda d : d.eqsys.T_cold[1:,0], *args)
    _, Ti = combine(lambda d : d.eqsys.W_i.getTemperature('D')[1:,0], *args)
    _, Ip = combine(lambda d : d.eqsys.I_p[1:,0], *args)
    _, Iw = combine(lambda d : d.eqsys.I_wall[1:,0], *args)
    a     = so1.grid.r_f[-1]
    Bphi  = 2.3
    Beddy = scipy.constants.mu_0 / (1.0*np.pi) * Iw
    Iref  = 1e5
    Bz    = 1e-3 + Beddy

    mD = scipy.constants.m_p + scipy.constants.m_n

    Lf = 0.25 * a * Bphi / Bz * np.exp(Ip / Iref)
    Cs = np.sqrt(scipy.constants.e*(Te+Ti)/mD)
    tauPar = Lf / Cs

    idx = 10000 + 237
    print(f'Te = {Te[idx]}')
    print(f'Ti = {Ti[idx]}')
    print(f'Ip = {Ip[idx]}')
    print(f'Iw = {Iw[idx]}')

    return t, tauPar


def getTauPerp(*args):
    so1 = args[0]

    t, Te = combine(lambda d : d.eqsys.T_cold[1:,0], *args)
    a = so1.grid.r_f[-1]
    Bphi  = 2.3

    DBohm = 1/16 * Te/Bphi
    vBohm = 2*DBohm / a

    return t, a / vBohm


so1 = STREAMOutput('output11.h5')
so2 = STREAMOutput('output12.h5')

_, Tcold_ion = combine(lambda d : d.other.fluid.Tcold_ion_coll[:,0], so1, so2)

"""
t, P_equi = getPequi(so1, so2)
plt.plot(t, Tcold_ion, 'k', label=r'Tcold\_ion\_coll')
plt.plot(t, P_equi, 'r--', label=r'$P_{\rm equi}$')
#plt.plot(t, P_equi / Tcold_ion, 'k')
plt.xlabel(r'Time (s)')
plt.ylabel('W/m$^-3$')
plt.legend()
"""

t, tauPar = getTauPar(so1, so2)
t, tauPerp = getTauPerp(so1, so2)

_, cTauPar  = combine(lambda d : d.other.stream.tau_D_par[:,0], so1, so2)
_, cTauPerp = combine(lambda d : d.other.stream.tau_D_perp[:,0], so1, so2)

_, axs = plt.subplots(1,2)

axs[0].plot(t, tauPar, 'k')
axs[0].plot(t, cTauPar, 'r--')
axs[1].plot(t, tauPerp, 'k')
axs[1].plot(t, cTauPerp, 'r--')

plt.show()

