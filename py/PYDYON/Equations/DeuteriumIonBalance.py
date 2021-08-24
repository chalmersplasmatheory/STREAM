# Particle balance equation for ionized deuterium 

from .. ADAS import ADAS


class DeuteriumIonBalance:
    

    def __init__(self, quantities, ions):
        """
        Constructor.
        """
        self.adas = ADAS()
        self.quantities = quantities
        self.ions = ions


    def __call__(self, t, x):
        return self.eval(t, x)


    def eval(self, t, x, full=False):
        """
        Evaluate the deuterium ion balance term.
        """
        Vn = self.quantities.getV_n('D')
        Vp = self.quantities.getV_p()

        ne = self.quantities['ne']
        Te = self.quantities['Te']
        Ti = self.quantities['Ti']
        nD = self.quantities['niD']

        Rrec = self.adas.ACD('D', 1, n=ne, T=Te)
        Riz  = self.adas.SCD('D', 0, n=ne, T=Te)

        # Ionization & recombination
        posIoniz = Vn/Vp * Riz*ne*nD[0]
        negIoniz = 0
        posRec = 0
        negRec = -Rrec*ne*nD[1]
        izrec = posIoniz + negRec

        # Charge exchange
        cx = 0
        for ion in self.ions:
            A = ion['name']
            if A == 'D':
                continue
            
            Z = ion['Z']
            ni = self.quantities.getIonData(A)

            for Z0 in range(1,Z+1):
                Rcx = self.adas.CCD(A, Z0, n=ni[Z0], T=Te)
                cx += Rcx*nD[0] * ni[Z0]

        dn = izrec + Vn/Vp * cx

        posCX = Vn/Vp*cx
        negCX = 0

        if full:
            return dn, posIoniz, negIoniz, posRec, negRec, posCX, negCX
        else:
            return dn


