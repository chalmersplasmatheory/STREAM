# Energy redistribution due to equilibration between electrons and ions

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


