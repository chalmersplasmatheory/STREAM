# Convective power loss term for electrons

import scipy.constants
from .. ConfinementTime import ConfinementTime


class ElectronConvectivePowerTerm:
    

    def __init__(self, quantities, ions, Bphi=2.3, Bv=1e-3, l_MK2=1):
        """
        Constructor.
        """
        self.quantities = quantities
        self.ions = ions

        self.tau = ConfinementTime(quantities, ions, Bphi=Bphi, Bv=Bv, l_MK2=l_MK2)


    def __call__(self, t, x):
        return self.eval(t, x)


    def eval(self, t, x):
        """
        Evaluate this ElectronConvectivePowerTerm.
        """
        ne = self.quantities['ne']
        Te = self.quantities['Te']
        e  = scipy.constants.e

        return 3/2 * ne*e*Te / self.tau(t, x)

