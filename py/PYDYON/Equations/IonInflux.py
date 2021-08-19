# Ion influx due to sputtering and recycling


from .. ConfinementTime import ConfinementTime


class IonInflux:
    

    def __init__(self, quantities, ions, Bphi=2.3, Bv=1e-3, l_MK2=1):
        """
        Constructor.

        :param quantities: UnknownQuantityHandler object.
        :param ions:       IonHandler object.
        """
        self.quantities = quantities
        self.ions = ions
        self.tau = ConfinementTime(quantities, ions, Bphi=Bphi, Bv=Bv, l_MK2=l_MK2)


    def __call__(self, x, ionname):
        return self.eval(x, ionname)


    def eval(self, x, ionname):
        """
        Evaluate this deuterium influx term.
        """
        Vp = self.quantities.getV_p()
        Vn_tot = self.quantities.getV_n_tot('D')

        Gamma = 0
        for ion in self.ions:
            A = ion['name']
            Z = ion['Z']

            ni = self.quantities.getIonData(A)
            for Z0 in range(1, Z+1):
                Y = self.evaluateRecyclingCoefficient(ionname, A)
                Gamma += Y*ni[Z0]

        Gamma *= Vp/(Vn_tot*self.tau(x))
        return Gamma


    def evaluateRecyclingCoefficient(I, A):
        """
        Evaluate the sputtering/recycling coefficient.
        """
        if I == 'D':
            if A == 'D': raise Exception('This term should not be applied to deuterium.')
            elif A == 'C': return 0
            elif A == 'O': return 0
        elif I == 'C':
            if A == 'D': return 0.03
            elif A == 'C': return 0
            elif A == 'O': return 1
        elif I == 'O':
            if A == 'D': return 0
            elif A == 'C': return 0
            elif A == 'O': return 1
        else:
            raise Exception(f'Unrecognized ion species: {I}.')

        raise Exception(f'Unrecognized ion species: {A}.')


