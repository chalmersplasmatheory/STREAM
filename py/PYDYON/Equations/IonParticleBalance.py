# Particle balance equation for ionized deuterium 

from .. ADAS import ADAS


class IonParticleBalance:
    

    def __init__(self, quantities, ions):
        """
        Constructor.
        """
        self.adas = ADAS()
        self.quantities = quantities
        self.ions = ions


    def __call__(self, x):
        return self.eval(x)


    def eval(self, x, ionname, ionZ0):
        """
        Evaluate the deuterium ion balance term.
        """
        Vn_tot = self.quantities.getV_n_tot(ionname)
        Vn = self.quantities.getV_n(ionname)
        VnD = self.quantities.getV_n(ionname)
        Vp = self.quantities.getV_p()

        ne = self.quantities['ne']
        Te = self.quantities['Te']

        dndt = 0
        if ionZ0 == 0:  # Neutral
            Rrec = self.adas.ACD(ionname, 1, n=ne, T=Te)
            Riz  = self.adas.SCD(ionname, 0, n=ne, T=Te)

            ni = self.quantities.getIonData(ionname)

            dndt = -Vn/Vn_tot * Riz * ne*ni[0] + Vp/Vn_tot * Rrec * ne*ni[1]
        else:    # Ionized
            nD = self.quantities.getIonData('D')

            Riz12  = self.adas.SCD(ionname, ionZ0,   n=ne, T=Te)
            Riz01  = self.adas.SCD(ionname, ionZ0-1, n=ne, T=Te)
            Rrec10 = self.adas.ACD(ionname, ionZ0,   n=ne, T=Te)
            Rrec21 = self.adas.ACD(ionname, ionZ0+1, n=ne, T=Te)
            Rcx10  = self.adas.CCD(ionname, ionZ0,   n=ne, T=Te)
            Rcx21  = self.adas.CCD(ionname, ionZ0+1, n=ne, T=Te)

            Vfac = Vn/Vp if ionZ0==1 else 1.0

            # Ionization
            dndt  = Vfac*Riz01*ne*ni[ionZ0-1] - Riz12*ne*ni[ionZ0]
            # Recombination
            dndt += Rrec21*ne*ni[ionZ0+1] - Rrec10*ne*ni[ionZ0]
            # Charge-exchange
            dndt += VnD/Vp * (Rcx21 * nD[0]*ni[ionZ0+1] - Rcx10 * nD[0]*ni[ionZ0])

        return dndt


