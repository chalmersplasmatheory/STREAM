#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys
from ITERbaseline import generate_baseline
from STREAM import STREAMSettings, runiface


def advance(ss, tMax, nt=1000, output='output.h5', chain=True):
    """
    Advance the given simulation to the specified 'tMax' with the
    specified number of time steps.
    """
    ss2 = STREAMSettings(ss, chain=chain)

    ss2.timestep.setTmax(tMax)
    ss2.timestep.setNt(nt)

    return ss2, runiface(ss2, output)


def DTWithRunaways(runaways=True, ext=''):
    """
    D-T simulation including runaways.
    """
    ss = generate_baseline(prefill=1e-6, tritium=True, runaways=runaways)
    pfx = 'w' if runaways else 'no'

    ss, so1 = advance(ss, tMax=1e-4, nt=1e4, output=f'output/DT{pfx}RE-{ext}-1.h5', chain=False)
    ss, so2 = advance(ss, tMax=8-ss.timestep.tmax, nt=1e5, output=f'output/DT{pfx}RE-{ext}-2.h5')


def DTWithRunawaysHighDensity():
    ss = generate_baseline(prefill=1.5e-5, tritium=True, runaways=True)

    pfx = 'w' #if runaways else 'no'

    ss, so1 = advance(ss, tMax=5e-3, nt=2e4, output=f'output/DT{pfx}RE-1.h5', chain=False)
    ss, so2 = advance(ss, tMax=8-ss.timestep.tmax, nt=3e5, output=f'output/DT{pfx}RE-2.h5')


def DTNoRunaways():
    """
    D-T simulation including runaways.
    """
    DTWithRunaways(False)


def DWithFueling(t0=0.5):
    ss = generate_baseline(prefill=1e-6, tritium=False, runaways=True)
    pfx = 'w'

    ext = 'fuel{:03d}'.format(int(t0*100))

    ss, so1 = advance(ss, tMax=1e-4, nt=1e4, output=f'output/D{pfx}RE-{ext}-1.h5', chain=False)

    t = np.linspace(0, 8-ss.timestep.tmax, 100000)
    f = np.zeros(t.shape)
    f[np.where((t >= t0) & (t <= (t0+2)))] = 5e17
    ss.eqsys.n_i.setFueling('D', f, times=t)

    ss, so2 = advance(ss, tMax=8-ss.timestep.tmax, nt=1e5, output=f'output/D{pfx}RE-{ext}-2.h5')


def main():
    #DTWithRunaways(ext='fuel300')
    #DTNoRunaways()
    #DTWithRunawaysHighDensity()

    DWithFueling(t0=0.5)
    DWithFueling(t0=1.0)
    DWithFueling(t0=3.0)

    return 0


if __name__ == '__main__':
    sys.exit(main())


