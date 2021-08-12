# Compare circuit equation solutions for DYON and STREAM parameters
# using loop voltage and plasma resistance from a STREAM simulation.

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../py')
from STREAM import STREAMOutput
import solver
import scipy.constants
from scipy.interpolate import interp1d


def run(b, tmax, nt=100):
    mu0 = scipy.constants.mu_0
    R = 3
    a = 0.9

    L = mu0*R*np.log(R/a)
    Lext = mu0*R*np.log(R/b)
    Mew = mu0*R*np.log(b/a)

    print(f'L = {L}, Lext = {Lext}, Mew = {Mew}')

    STREAM = {'M': Lext, 'L_MK2': Lext, 'L_p': (L+Lext+Mew), 'R_MK2': 3.4e-4}

    sol = solver.solve(**STREAM, V=V, R_p=R_p, tmax=tmax, nt=nt)
    
    return sol.y[0][-1]


def scanb(V, R_p, tmax, nt=100):
    b = np.linspace(0.2, 1.0, 20)
    Ip = []
    for _b in b:
        Ip.append(run(_b, tmax=tmax, nt=nt))

    plt.plot(b, Ip)
    plt.xlabel('$b$')
    plt.ylabel(r'$I_{\rm p}$')
    plt.show()


so = STREAMOutput('../DYON_tests/output_initial.h5')

Rm = 3  # Tokamak major radius
a  = 0.9  # Tokamak minor radius
R0 = 2*np.pi*Rm/(np.pi*a**2)

sigma = np.insert(so.other.fluid.conductivity[:,0], 0, so.other.fluid.conductivity[0,0])

V = interp1d(so.grid.t[:], so.eqsys.V_loop_trans[:,0])
R_p = interp1d(so.grid.t[:], R0/sigma)

solver.compare(V=V, R_p=R_p, tmax=so.grid.t[-1])
#scanb(V=V, R_p=R_p, tmax=so.grid.t[-1])

