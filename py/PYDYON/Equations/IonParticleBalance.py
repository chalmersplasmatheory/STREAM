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


    def __call__(self, t, x, ionname, Z0):
        return self.eval(t, x, ionname, Z0)


    def eval(self, t, x, ionname, Z0):
        """
        Evaluate the deuterium ion balance term.
        """
        Vn_tot = self.quantities.getV_n_tot(ionname)
        Vn = self.quantities.getV_n(ionname)
        VnD = self.quantities.getV_n(ionname)
        Vp = self.quantities.getV_p()

        ne = self.quantities['ne']
        Te = self.quantities['Te']
        Ti = self.quantities['Ti']

        dndt = 0
        if Z0 == 0:  # Neutral
            Rrec = self.adas.ACD(ionname, 1, n=ne, T=Te)
            Riz  = self.adas.SCD(ionname, 0, n=ne, T=Te)

            ni = self.quantities.getIonData(ionname)

            dndt = -Vn/Vn_tot * Riz * ne*ni[0] + Vp/Vn_tot * Rrec * ne*ni[1]
        else:    # Ionized
            nD = self.quantities.getIonData('D')

            Riz12  = self.adas.SCD(ionname, Z0,   n=ne, T=Te)
            Riz01  = self.adas.SCD(ionname, Z0-1, n=ne, T=Te)
            Rrec10 = self.adas.ACD(ionname, Z0,   n=ne, T=Te)
            Rrec21 = self.adas.ACD(ionname, Z0+1, n=ne, T=Te)
            Rcx10  = self.adas.CCD(ionname, Z0,   n=ni[Z0], T=Ti)
            Rcx21  = self.adas.CCD(ionname, Z0+1, n=ni[Z0], T=Ti)

            Vfac = Vn/Vp if Z0==1 else 1.0

            # Ionization
            dndt  = Vfac*Riz01*ne*ni[Z0-1] - Riz12*ne*ni[Z0]
            # Recombination
            dndt += Rrec21*ne*ni[Z0+1] - Rrec10*ne*ni[Z0]
            # Charge-exchange
            dndt += VnD/Vp * (Rcx21 * nD[0]*ni[Z0+1] - Rcx10 * nD[0]*ni[Z0])

        return dndt


