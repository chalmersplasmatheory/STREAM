# Various helper routines for evaluating plasma and neutral volumes

from . ADAS import ADAS
import numpy as np


class PlasmaVolume:

    def __init__(self, a, R, V_vessel, ions, kappa=1):
        """
        Constructor.

        :param a:        Plasma minor radius.
        :param R:        Plasma major radius.
        :param V_vessel: Total vacuum vessel volume.
        :param ions:     IonHandler object.
        :param kappa:    Plasma elongation.
        """
        self.adas = ADAS()
        self.a = a
        self.R = R
        self.V_vessel = V_vessel
        self.ions = ions
        self.kappa = kappa


    def getV_p(self):
        """
        Evaluate the plasma volume V_p.
        """
        return 2*np.pi**2 * self.R * self.kappa*self.a**2


    def getLambda(self, ion, ne, Te, Ti):
        """
        Calculates the mean-free path length lambda_i for ion
        species 'ion'.

        :param ion: Name of ion species to calculate mean-free path for.
        :param ne:  Electron density to evaluate at.
        :param Te:  Electron temperature to evaluate at.
        :param Ti:  Ion temperature to evaluate at.
        """
        m = self.ions.getMass(ion)

        return np.sqrt(2*Ti/m) / (ne*self.adas.SCD(ion, 0, n=ne, T=Te))


    def getV_n(self, ion, ne, Te, Ti):
        """
        Calculates the plasma neutral volume occupied by ions of
        species 'ion'.

        :param ion: Name of ion species to calculate mean-free path for.
        :param ne:  Electron density to evaluate at.
        :param Te:  Electron temperature to evaluate at.
        :param Ti:  Ion temperature to evaluate at.
        """
        lambdai = self.getLambda(ion, ne=ne, Te=Te, Ti=Ti)

        if lambdai > self.a:
            return self.getV_p()
        else:
            return 2*np.pi*R*(np.pi*self.kappa*self.a**2 - np.pi*self.kappa*(a-lambdai)**2)


    def getV_n_tot(self, ion, ne, Te, Ti):
        """
        Calculates the *total* volume occupied by neutrals
        of the named ion species.

        :param ion: Name of ion species to calculate mean-free path for.
        :param ne:  Electron density to evaluate at.
        :param Te:  Electron temperature to evaluate at.
        :param Ti:  Ion temperature to evaluate at.
        """
        Vp = self.getV_p()
        Vn = self.getV_n(ion, ne=ne, Te=Te, Ti=Ti)

        return self.V_vessel - (Vp-Vn)


