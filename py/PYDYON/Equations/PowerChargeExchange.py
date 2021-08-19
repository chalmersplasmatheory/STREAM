# Charge exchange power loss term


from .. ADAS import ADAS


class ChargeExchangePowerTerm:
    

    def __init__(self, quantities, ions):
        """
        Constructor.
        """
        self.quantities = quantities
        self.ions = ions
        self.adas = ADAS()


    def __call__(self, x):
        return self.eval(x)


    def eval(self, x):
        """
        Evaluate the power lost through charge-exchange.
        """
        Vn = self.quantities.getV_n('D')
        Vp = self.quantities.getV_p()
        e  = scipy.constants.e

        nD0 = self.quantities['nD'][0]
        Ti  = self.quantities['Ti']
        ne  = self.quantities['ne']
        Te  = self.quantities['Te']
        T0  = 0.026

        rad = 0
        for ion in ions:
            A = ion['name']
            ni1 = self.quantities.getIonData(A)[1]

            rad += self.adas.CCD(A, 0, n=ne, T=Te) * ni1

        Pcx = Vn/Vp * (3/2) * nD0*e*(Ti-T0) * rad

        return Pcx
            

