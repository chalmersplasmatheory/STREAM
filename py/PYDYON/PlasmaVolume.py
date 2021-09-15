# Various helper routines for evaluating plasma and neutral volumes

from . ADAS import ADAS
import numpy as np
import scipy.interpolate
import scipy.constants


class PlasmaVolume:

    def __init__(self, a, R, V_vessel, ions, kappa=1, ta=None):
        """
        Constructor.

        :param a:        Plasma minor radius.
        :param R:        Plasma major radius.
        :param V_vessel: Total vacuum vessel volume.
        :param ions:     IonHandler object.
        :param kappa:    Plasma elongation.
        """
        self.adas = ADAS()

        if isinstance(a,np.ndarray):
            self.a = scipy.interpolate.interp1d(ta, a, bounds_error=False, fill_value="extrapolate")
        elif not callable(a):
            self.a = lambda ta : a
        else:
            self.a = a

        self.R = R
        self.V_vessel = V_vessel
        self.ions = ions
        self.kappa = kappa


    def getV_p(self, t):
        """
        Evaluate the plasma volume V_p.
        """
        return 2*np.pi**2 * self.R * self.kappa*self.a(t)**2


    def getLambda(self, ion, ne, Te, Ti):
        """
        Calculates the mean-free path length lambda_i for ion
        species 'ion'.

        :param ion: Name of ion species to calculate mean-free path for.
        :param ne:  Electron density to evaluate at.
        :param Te:  Electron temperature to evaluate at.
        :param Ti:  Ion temperature to evaluate at.
        """
        e = scipy.constants.e
        m = self.ions.getMass(ion)

        return np.sqrt(2*e*Ti/m) / (ne*self.adas.SCD(ion, 0, n=ne, T=Te))


    def getV_n(self, t, ion, ne, Te, Ti):
        """
        Calculates the plasma neutral volume occupied by ions of
        species 'ion'.

        :param ion: Name of ion species to calculate mean-free path for.
        :param ne:  Electron density to evaluate at.
        :param Te:  Electron temperature to evaluate at.
        :param Ti:  Ion temperature to evaluate at.
        """
        lambdai = self.getLambda(ion, ne=ne, Te=Te, Ti=Ti)

        if lambdai > self.a(t):
            return self.getV_p(t)
        else:
            R = self.R
            a = self.a(t)
            pi = np.pi
            kappa = self.kappa
            return 2*pi*R*(pi*kappa*a**2 - pi*kappa*(a-lambdai)**2)


    def getV_n_tot(self, t, ion, ne, Te, Ti):
        """
        Calculates the *total* volume occupied by neutrals
        of the named ion species.

        :param ion: Name of ion species to calculate mean-free path for.
        :param ne:  Electron density to evaluate at.
        :param Te:  Electron temperature to evaluate at.
        :param Ti:  Ion temperature to evaluate at.
        """
        Vp = self.getV_p(t)
        Vn = self.getV_n(t, ion, ne=ne, Te=Te, Ti=Ti)

        return self.V_vessel - (Vp-Vn)


