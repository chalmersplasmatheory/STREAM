# Transport of ions

from .. ConfinementTime import ConfinementTime


class IonTransport:
    

    def __init__(self, quantities, ions, Bphi=2.3, Bv=1e-3, l_MK2=1):
        """
        Constructor.

        :param quantities: UnknownQuantityHandler object.
        :param ions:       IonHandler object.
        """
        self.tau = ConfinementTime(quantities, ions, Bphi=Bphi, Bv=Bv, l_MK2=l_MK2)
        self.ions = ions
        self.quantities = quantities


    def __call__(self, t, x, ionname, Z0):
        return self.eval(t, x, ionname, Z0)


    def eval(self, t, x, ionname, Z0):
        """
        Evaluate this ion transport term.

        :param x:       Latest solution of equations.
        :param ionname: Name of ion to apply transport to.
        :param Z0:      Charge state of ion to apply transport to.
        """
        niZ = self.quantities.getIonData(ionname)[Z0]

        return niZ / self.tau(t, x)


