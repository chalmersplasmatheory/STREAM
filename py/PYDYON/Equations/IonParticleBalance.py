# Particle balance equation for ionized deuterium 

from .. ADAS import ADAS
import numpy as np
import sys


class IonParticleBalance:
    

    def __init__(self, quantities, ions):
        """
        Constructor.
        """
        self.adas = ADAS()
        self.quantities = quantities
        self.ions = ions
        self.mainIon = ions.getMainIonName()


    def __call__(self, t, x, ionname, Z0):
        return self.eval(t, x, ionname, Z0)


    def eval(self, t, x, ionname, Z0):
        """
        Evaluate the deuterium ion balance term.
        """
        Vn_tot = self.quantities.getV_n_tot(ionname)
        Vn = self.quantities.getV_n(ionname)
        VnD = self.quantities.getV_n(self.mainIon)
        Vp = self.quantities.getV_p()

        ne = self.quantities['ne']
        Te = self.quantities['Te']
        Ti = self.quantities['Ti']
        nD = self.quantities.getIonData(self.mainIon)

        dndt = 0
        if Z0 == 0:  # Neutral
            #print(f'A={ionname}-{Z0}: ne={ne}, Te={Te}')
            Rrec = self.adas.ACD(ionname, 1, n=ne, T=Te)
            Riz  = self.adas.SCD(ionname, 0, n=ne, T=Te)
            #print(f'    Rrec={Rrec}, Riz={Riz}')

            ni = self.quantities.getIonData(ionname)

            Rcx = self.adas.CCD(ionname, 1, n=ni[1], T=Ti)

            dndt = -Vn/Vn_tot * Riz * ne*ni[0] + Vp/Vn_tot * Rrec * ne*ni[1] + VnD/Vn_tot * Rcx * nD[0] * ni[1]
        else:    # Ionized
            ni = self.quantities.getIonData(ionname)

            if Z0 >= self.ions[ionname]['Z']:
                Riz12 = 0
                Rrec21 = 0
                Rcx21 = 0
                nip1 = 0
            else:
                Riz12  = self.adas.SCD(ionname, Z0,   n=ne, T=Te)
                Rrec21 = self.adas.ACD(ionname, Z0+1, n=ne, T=Te)
                Rcx21  = self.adas.CCD(ionname, Z0+1, n=ni[Z0], T=Ti)
                nip1 = ni[Z0+1]

            Riz01  = self.adas.SCD(ionname, Z0-1, n=ne, T=Te)
            Rrec10 = self.adas.ACD(ionname, Z0,   n=ne, T=Te)
            Rcx10  = self.adas.CCD(ionname, Z0,   n=ni[Z0], T=Ti)

            Vfac = Vn/Vp if Z0==1 else 1.0

            # Ionization
            dndt  = Vfac*Riz01*ne*ni[Z0-1] - Riz12*ne*ni[Z0]
            # Recombination
            dndt += Rrec21*ne*nip1 - Rrec10*ne*ni[Z0]
            # Charge-exchange
            dndt += VnD/Vp * (Rcx21 * nD[0]*nip1 - Rcx10 * nD[0]*ni[Z0])

        return dndt


