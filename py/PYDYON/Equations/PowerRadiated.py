# Radiated + ionization power term


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


    def __call__(self, t, x):
        return self.eval(t, x)


    def eval(self, t, x):
        """
        Evaluate the radiated power.
        """
        e = scipy.constants.e
        ne = self.quantities['ne']
        Te = self.quantities['Te']
        Vp = self.quantities.getPlasmaVolume()
        
        P = 0
        for ion in self.ions:
            Z = ion['Z']
            A = ion['name']

            ni = self.quantities.getIonData(ion)
            Vn = self.quantities.getNeutralVolume(A)

            # Z0=0 contribution
            P += Vn/Vp * self.adas.PLT(A, 0, n=ne, T=Te) * ne * ni[0]
            # R_iz
            P += Vn/Vp * self.adas.SCD(A, 0, n=ne, T=Te) * e*self.nist.ionization(A, 0) * ne * ni[0]
            
            # Z0>0 contributions
            for Z0 in range(1, Z+1):
                q = 0
                if Z0 < Z:
                    # P_line
                    q += self.adas.PLT(A, Z0, n=ne, T=Te)
                    # R_iz
                    q += self.adas.SCD(A, Z0, n=ne, T=Te) * self.nist.ionization(A, Z0) * e

                # P_RB
                q += self.adas.PRB(A, Z0, n=ne, T=Te)
                # R_rec
                q -= self.adas.ACD(A, Z0, n=ne, T=Te) * self.nist.ionization(A, Z0-1) * e

                P += q*ne*ni[Z0]

        return P


