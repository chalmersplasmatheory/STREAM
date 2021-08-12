# Routines for solving the DYON circuit equations

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
from scipy.integrate import solve_ivp


def solve(M, L_MK2, L_p, R_MK2, V, R_p, tmax=0.5, nt=100, Ip0=0):
    """
    Solve the DYON circuit equations.

    :param M:     Mutual inductance between plasma and MK2 ring.
    :param L_MK2: Self-inductance of MK2 ring.
    :param L_p:   Self-inductance of plasma.
    :param R_MK2: Resistance of MK2 ring.
    :param V:     Applied loop voltage.
    :param R_p:   Plasma resistance.
    """
    _V = V
    _R_p = R_p

    if not callable(V): _V = lambda t : V
    if not callable(R_p): _R_p = lambda t : R_p

    def dIdt(t, I):
        Ip = I[0]
        IMK2 = I[1]

        pf = M - L_MK2*L_p/M
        eq1 = (_V(t) - R_MK2 * IMK2 - L_MK2/M * (_V(t) - _R_p(t)*Ip)) / pf
        eq2 = (_V(t) - _R_p(t) * Ip  - L_p/M *  (_V(t) - R_MK2*IMK2)) / pf

        return [eq1, eq2]

    return solve_ivp(dIdt, [0 , tmax], [Ip0, 0], t_eval=np.linspace(0, tmax, nt))


def compare(V, R_p, tmax=0.5, nt=100):
    DYON   = {'M': 2.49e-6, 'L_MK2': 9.1e-6, 'L_p': 5.4e-6, 'R_MK2': 7.5e-4}
    STREAM = {'M': 4.17e-6, 'L_MK2': 4.17e-6, 'L_p': 8.3e-6, 'R_MK2': 3.4e-4}
    #STREAM = {'M': 4.17e-6, 'L_MK2': 9.1e-6, 'L_p': 8.3e-6, 'R_MK2': 3.4e-4}   # better

    sol_DYON = solve(**DYON, V=V, R_p=R_p, tmax=tmax, nt=nt)
    sol_STREAM = solve(**STREAM, V=V, R_p=R_p, tmax=tmax, nt=nt)

    plt.plot(sol_DYON.t, sol_DYON.y[0], 'r-', label=r'$I_{\rm p}$')
    plt.plot(sol_DYON.t, sol_DYON.y[1], 'k-', label=r'$I_{\rm MK2}$')
    plt.plot(sol_STREAM.t, sol_STREAM.y[0], 'r--', label='STREAM')
    plt.plot(sol_STREAM.t, sol_STREAM.y[1], 'k--')

    Imax = max(np.amax(sol_DYON.y), np.amax(sol_STREAM.y))
    plt.axis([0, np.amax(sol_DYON.t), 0, 1.1*Imax])
    plt.ylabel('Current (A)')

    plt.legend()
    plt.show()


if __name__ == '__main__':
    V = 11
    R_p = 1e-5

    compare(V, R_p)


