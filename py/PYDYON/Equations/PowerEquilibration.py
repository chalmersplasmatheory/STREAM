# Energy redistribution due to equilibration between electrons and ions

import numpy as np
import scipy.constants
from .. import Conductivity


class EquilibrationPowerTerm:
    

    def __init__(self, quantities, ions):
        """
        Constructor.

        :param quantities: UnknownQuantityHandler object.
        :param ions:       IonHandler object.
        """
        self.quantities = quantities
        self.ions = ions


    def __call__(self, t, x):
        return self.eval(t, x)


    def eval(self, t, x):
        """
        Evaluate the equilibration power.
        """
        #return self.evalDYON(t, x)
        return self.evalSTREAM(t,x)


    def evalDYON(self, t, x):
        """
        Evaluate this term as it is implemented in DYON.
        """
        ne = self.quantities['ne']
        Te = self.quantities['Te']
        Ti = self.quantities['Ti']
        logLambda = Conductivity.getCoulombLogarithm(T=Te, n=ne)
        #logLambda = 10
        amu = scipy.constants.physical_constants['atomic mass constant'][0]

        pf = 7.75e-34 * (Te - Ti) * (ne*logLambda)/Te**1.5

        nz = 0
        for ion in self.ions:
            Z  = ion['Z']
            A  = ion['name']
            ni = self.quantities.getIonData(A)

            for Z0 in range(1,Z+1):
                nz += ni[Z0] * Z0**2 / (self.ions.getMass(A) / amu)

        return pf*nz


    def evalSTREAM(self, t, x):
        """
        Evaluate this term as it is implemented in STREAM.
        """
        ne = self.quantities['ne']
        Te = self.quantities['Te']
        Ti = self.quantities['Ti']
        logLambda = Conductivity.getCoulombLogarithm(T=Te, n=ne)

        e = scipy.constants.e
        eps0 = scipy.constants.epsilon_0
        me = scipy.constants.m_e
        mi = self.ions.getMass(self.ions.getMainIonName())

        pf = e**4 * ne * logLambda / (((2*np.pi)**1.5)*(eps0**2)*me) / np.sqrt(e)

        nz = 0
        for ion in self.ions:
            Z  = ion['Z']
            A  = ion['name']
            ni = self.quantities.getIonData(A)

            if Z != 1:
                continue

            for Z0 in range(1,Z+1):
                nz += ni[Z0] * Z0**2 / (self.ions.getMass(A))
        
        Tfac = (Te-Ti)/(Te/me + Ti/mi)**1.5

        return pf*nz*Tfac


