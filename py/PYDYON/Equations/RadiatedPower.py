# Radiated power term


import scipy.constants
from .. ADAS import ADAS
from .. NIST import NIST


class RadiatedPowerTerm:

    
    def __init__(self, quantities, ions):
        """
        Constructor.

        :param quantities: Helper class for extracting useful quantities.
        :param ions:       Definition of ions included in simulation.
        """
        self.adas = ADAS() 
        self.nist = NIST()
        self.ions = ions
        self.quantities = quantities


    def Prad(self, x):
        """
        Evaluate the radiated power.
        """
        ne = self.quantities['ne']
        Te = self.quantities['Te']
        Vp = self.quantities.getPlasmaVolume()
        
        P = 0
        for ion in self.ions:
            Z = ion['Z']
            A = ion['name']

            ni = self.quantities.getIonData(ion)
            Vn = self.quantities.getNeutralVolume(A)

            P += Vn/Vp * self.adas.PLT(A, 0, n=ne, T=Te) * ne * ni[0]
            
            for Z0 in range(1, Z+1):
                Wioniz = self.nist.ionization(A, Z0-1) * scipy.constants.e

                q = 0
                if Z0 < Z:
                    q  = self.adas.PLT(A, Z0, n=ne, T=Te)

                q += self.adas.PRB(A, Z0, n=ne, T=Te)
                q -= self.adas.ACD(A, Z0, n=ne, T=Te) * Wioniz

                P += q*ne*ni[Z0]

        return P


