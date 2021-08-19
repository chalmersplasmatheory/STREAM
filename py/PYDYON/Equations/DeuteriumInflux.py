# Deuterium influx term


class DeuteriumInflux:


    # Parameters of deuterium recycling coefficient
    c1 = 1.1
    c2 = 0.09
    c3 = 0.1
    

    def __init__(self, quantities, ions, simple=False, Bphi=2.3, Bv=1e-3, l_MK2=1):
        """
        Constructor.

        :param quantities: UnknownQuantityHandler object.
        :param ions:       IonHandler object.
        :param simple:     If ``True``, assumes that the recycling coefficient is Y_D(t) = 1.
        """
        self.quantities = quantities,
        self.ions = ions
        self.simpleYD = simple
        self.tau = ConfinementTime(quantities, ions, Bphi=Bphi, Bv=Bv, l_MK2=l_MK2)


    def __call__(self, x):
        return self.eval(x)


    def eval(self, t, x):
        """
        Evaluate this deuterium influx term.
        """
        YD = self.evaluateRecyclingCoefficient(t)
        Vp = self.quantities.getV_p()
        Vn_tot = self.quantities.getV_n_tot('D')

        nD1 = self.quantities.getIonData('D')

        Gamma = Vp*YD*nD1 / self.tau(x)

        return Gamma / Vn_tot


    def evaluateRecyclingCoefficient(t):
        """
        Evaluate the deuterium recycling coefficient.
        """
        if self.simpleYD:
            return 1.0

        YD = self.c1 - self.c2*(1 - np.exp(-t/self.c3))
        return YD


