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
        self.mainIon = ions.getMainIonName()


    def __call__(self, t, x, ionname):
        return self.eval(t, x, ionname)


    def eval(self, t, x, ionname):
        """
        Evaluate this deuterium influx term.
        """
        Vp = self.quantities.getV_p()
        Vn_tot = self.quantities.getV_n_tot(ionname)

        Gamma = 0
        for ion in self.ions:
            A = ion['name']
            Z = ion['Z']

            ni = self.quantities.getIonData(A)
            for Z0 in range(1, Z+1):
                Y = self.evaluateRecyclingCoefficient(ionname, A)
                Gamma += Y*ni[Z0]

        #print(f'{ionname}:  Gamma = {Gamma}')
        Gamma *= Vp/(Vn_tot*self.tau(t, x))

        return Gamma


    def evaluateRecyclingCoefficient(self, I, A):
        """
        Evaluate the sputtering/recycling coefficient.
        """
        if I == self.mainIon:
            if A == self.mainIon: raise Exception('This term should not be applied to deuterium.')
            elif A == 'C': return 0
            elif A == 'O': return 0
        elif I == 'C':
            #if A == 'D': return 0.03
            if A == self.mainIon: return 0.015
            elif A == 'C': return 0
            elif A == 'O': return 1
        elif I == 'O':
            if A == self.mainIon: return 0
            elif A == 'C': return 0
            elif A == 'O': return 1
        else:
            raise Exception(f'Unrecognized ion species: {I}.')

        raise Exception(f'Unrecognized ion species: {A}.')


