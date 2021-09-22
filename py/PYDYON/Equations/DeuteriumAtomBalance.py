# Particle balance equation for neutral deuterium 

import numpy as np
from .. ADAS import ADAS


class DeuteriumAtomBalance:
    

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
        Evaluate the neutral deuterium particle balance term.

        :param full: If ``True``, also returns values of individual terms.
        """
        V_n_tot = self.quantities.getV_n_tot('D')
        Vn = self.quantities.getV_n('D')
        Vp = self.quantities.getV_p()

        ne = self.quantities['ne']
        Te = self.quantities['Te']
        Ti = self.quantities['Ti']
        nD = self.quantities['niD']

        Rrec = self.adas.ACD('D', 1, n=ne, T=Te)
        Riz  = self.adas.SCD('D', 0, n=ne, T=Te)

        # Ionization & recombination
        posRec   = Vp*Rrec*ne*nD[1]
        negRec   = 0
        posIoniz = 0
        negIoniz = -Vn*Riz*ne*nD[0]

        izrec = posRec + negIoniz

        # Charge exchange
        cx = 0
        for ion in self.ions:
            A = ion['name']
            if A == 'D':
                continue

            Z = ion['Z']
            ni = self.quantities.getIonData(A)

            for Z0 in range(1,Z+1):
                Rcx = self.adas.CCD(A, Z0, n=ni[Z0], T=Ti)
                cx += Rcx*nD[0] * ni[Z0]

                if np.isnan(cx):
                    print(f'Rcx  = {Rcx}')
                    print(f'Ti   = {Ti}')
                    print(f'niZ0 = {ni[Z0]}')

        negCX = -Vn*cx
        posCX = 0

        dn = 1/V_n_tot * (izrec + Vn * cx)

        if full:
            return dn, posIoniz/V_n_tot, negIoniz/V_n_tot, posRec/V_n_tot, negRec/V_n_tot, posCX/V_n_tot, negCX/V_n_tot
        else:
            return dn

    def evalTi(self,t,x):
        Ti = self.quantities['Ti']
        return Ti

    def evalniD(self,t,x):
        ni = self.quantities.getIonData('D')
        return ni[0], ni[1]

    def evalniC(self,t,x):
        ni = self.quantities.getIonData('C')
        return ni[0], ni[1], ni[2], ni[3], ni[4], ni[5], ni[6]

    def evalniO(self,t,x):
        ni = self.quantities.getIonData('O')
        return ni[0], ni[1], ni[2], ni[3], ni[4], ni[5], ni[6], ni[7], ni[8]
