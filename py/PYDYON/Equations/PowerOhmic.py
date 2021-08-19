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


    def __call__(self, x):
        return self.eval(x)


    def eval(self, x):
        """
        Evaluate the radiated power.
        """
        Ip = self.quantities['Ip']
        Vp = self.quantities.getPlasmaVolume()
        sg = Conductivity.eval(self.quantities)

        R = self.quantities.plasmavolume.R
        a = self.quantities.plasmavolume.a

        Rp = 2*R/a**2 * (1/sg)

        Pohm = Ip**2 * Rp / Vp
        return Pohm


