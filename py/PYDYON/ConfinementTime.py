# Evaluation of the particle confinement time

import numpy as np
import scipy.constants


class ConfinementTime:
    

    def __init__(self, quantities, ions, Bphi=2.3, Bv=1e-3, l_MK2=1):
        """
        Constructor.

        :param quantities: UnknownQuantityHandler object.
        :param ions:       IonHandler object.
        :param Bphi:       Toroidal magnetic field strength (T).
        :param Bv:         Verical magnetic field strength (T).
        :param l_MK2:      Distance between plasma centre and passive structure (m).
        """
        self.Bphi = Bphi
        self.Bv = Bv
        self.l_MK2 = l_MK2
        self.quantities = quantities
        self.ions = ions


    def __call__(self, t, x):
        return self.eval(t, x)


    def eval(self, t, x):
        """
        Evaluate this confinement time.
        """
        tauPar  = self.evalParallel(t, x)
        tauPerp = self.evalPerpendicular(t, x)
        itau = 1/tauPar + 1/tauPerp

        return 1/itau


    def evalParallel(self, t, x):
        """
        Evaluate parallel confinement time.
        """
        Ip = self.quantities['Ip']
        IMK2 = self.quantities['IMK2']
        Te = self.quantities['Te']
        Ti = self.quantities['Ti']

        a    = self.quantities.plasmavolume.a(t)
        Bphi = self.Bphi
        Iref = 100e3
        e    = scipy.constants.e

        Beddy = scipy.constants.mu_0 * IMK2 / (np.pi * self.l_MK2)
        Bz = self.Bv + Beddy

        Lf = 0.25*a*Bphi/Bz * np.exp(Ip/Iref)

        mD = scipy.constants.m_p + scipy.constants.m_n
        Cs = np.sqrt(e*(Te + Ti)/mD)

        return Lf / Cs


    def evalPerpendicular(self, t, x):
        """
        Evaluate perpendicular confinement time.
        """
        Te = self.quantities['Te']

        a = self.quantities.plasmavolume.a(t)
        DBohm = 1/16 * Te/self.Bphi
        vBohm = 2*DBohm / a

        return a / vBohm


