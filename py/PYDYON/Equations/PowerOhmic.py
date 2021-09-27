# Ohmic power term


import scipy.constants
from .. import Conductivity


class OhmicPowerTerm:

    
    def __init__(self, quantities, ions):
        """
        Constructor.

        :param quantities: UnknownQuantityHandler object.
        :param ions:       IonHandler object.
        """
        self.ions = ions
        self.quantities = quantities


    def __call__(self, t, x):
        return self.eval(t, x)


    def eval(self, t, x):
        """
        Evaluate the radiated power.
        """
        Ip = self.quantities['Ip']
        Vp = self.quantities.getV_p()
        Rp = Conductivity.evalResistance(self.quantities)

        Pohm = Ip**2 * Rp / Vp
        return Pohm


