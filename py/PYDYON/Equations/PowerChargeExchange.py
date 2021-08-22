# Charge exchange power loss term


import scipy.constants
from .. ADAS import ADAS


class ChargeExchangePowerTerm:
    

    def __init__(self, quantities, ions):
        """
        Constructor.
        """
        self.quantities = quantities
        self.ions = ions
        self.adas = ADAS()


    def __call__(self, t, x):
        return self.eval(t, x)


    def eval(self, t, x):
        """
        Evaluate the power lost through charge-exchange.
        """
        Vn = self.quantities.getV_n('D')
        Vp = self.quantities.getV_p()
        e  = scipy.constants.e

        nD0 = self.quantities['niD'][0]
        Ti  = self.quantities['Ti']
        ne  = self.quantities['ne']
        Te  = self.quantities['Te']
        T0  = 0.026

        rad = 0
        for ion in self.ions:
            A = ion['name']
            ni1 = self.quantities.getIonData(A)[1]

            #rad += self.adas.CCD(A, 1, n=ni1, T=Ti) * ni1
            rad += 1.066e-14*(Ti**.327) * ni1
            #print('ADAS: {},  Fit: {}'.format(self.adas.CCD(A, 1, n=ni1, T=Ti), 1e-14*Ti**.327))

        #Pcx = Vn/Vp * (3/2) * nD0*e*(Ti-T0) * rad
        Pcx = Vn/Vp * 1.5 * nD0*e*(Ti-T0) * rad

        return Pcx
            

