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
        self.mainIon = ions.getMainIonName()


    def __call__(self, t, x, full=False):
        return self.eval(t, x, full)


    def eval(self, t, x, full=False):
        """
        Evaluate the deuterium ion balance term.
        """
        Vn = self.quantities.getV_n(self.mainIon)
        Vp = self.quantities.getV_p()

        ne = self.quantities['ne']
        Te = self.quantities['Te']
        Ti = self.quantities['Ti']
        nD = self.quantities[f'ni{self.mainIon}']

        Rrec = self.adas.ACD(self.mainIon, 1, n=ne, T=Te)
        Riz  = self.adas.SCD(self.mainIon, 0, n=ne, T=Te)

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
            if A == self.mainIon:
                continue
            
            Z = ion['Z']
            ni = self.quantities.getIonData(A)

            for Z0 in range(1,Z+1):
                Rcx = self.adas.CCD(A, Z0, n=ni[Z0], T=Ti)
                cx += Rcx*nD[0] * ni[Z0]

        dn = izrec + Vn/Vp * cx

        posCX = Vn/Vp*cx
        negCX = 0

        if full:
            return dn, posIoniz, negIoniz, posRec, negRec, posCX, negCX
        else:
            return dn


