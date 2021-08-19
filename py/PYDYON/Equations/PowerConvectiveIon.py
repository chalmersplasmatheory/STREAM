# Convective power loss term for ions

import scipy.constants
from .. ConfinementTime import ConfinementTime


class IonConvectivePowerTerm:
    

    def __init__(self, quantities, ions, Bphi=2.3, Bv=1e-3, l_MK2=1):
        """
        Constructor.
        """
        self.quantities = quantities
        self.ions = ions

        self.tau = ConfinementTime(quantities, ions, Bphi=Bphi, Bv=Bv, l_MK2=l_MK2)


    def __call__(self, x):
        return self.eval(x)


    def eval(self, x):
        """
        Evaluate this IonConvectivePowerTerm.
        """
        Te = self.quantities['Ti']
        e  = scipy.constants.e

        N = 0
        for ion in self.ions:
            A = ion['name']
            Z = ion['Z']

            for Z0 in range(1, Z+1):
                ni = self.quantities.getIonData(A)

                P += ni[Z0]

        return 3/2 * N * e*Ti / self.tau(x)


