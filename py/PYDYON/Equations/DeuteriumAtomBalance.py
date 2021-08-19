# Particle balance equation for neutral deuterium 

from .. ADAS import ADAS


class DeuteriumAtomBalance:
    

    def __init__(self, quantities, ions):
        """
        Constructor.
        """
        self.adas = ADAS()
        self.quantities = quantities
        self.ions = ions


    def __call__(self, x):
        return self.eval(x)


    def eval(self, x):
        """
        Evaluate the neutral deuterium particle balance term.
        """
        V_n_tot = self.quantities.getV_n_tot('D')
        Vn = self.quantities.getV_n('D')
        Vp = self.quantities.getV_p()

        ne = self.quantities['ne']
        Te = self.quantities['Te']
        nD = self.quantities['niD']

        Rrec = self.adas.ACD('D', 1, n=ne, T=Te)
        Riz  = self.adas.SCD('D', 0, n=ne, T=Te)

        # Ionization & recombination
        izrec = Vp*Rrec*ne*nD[1] - Vn*Riz*ne*nD[0]

        # Charge exchange
        cx = 0
        for ion in self.ions:
            A = ion['name']
            if A == 'D':
                continue
            
            Z = ion['Z']
            ni = self.quantities.getIonData(A)

            for Z0 in range(1,Z+1):
                Rcx = self.adas.CCD(A, Z0, n=ne, T=Te)
                cx += Rcx*nD[0] * ni

        return 1/V_n_tot * (izrec - Vn * cx)


