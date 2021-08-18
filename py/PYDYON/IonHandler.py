# Ion handler

import scipy.constants


class IonHandler:
    

    def __init__(self):
        """
        Constructor.
        """
        self.ions = []


    def addIon(self, name, Z, m=None):
        """
        Add an ion to the simulation.

        :param name: Name of ion.
        :param Z:    Ion charge.
        :param n:    Ion initial densities.
        """
        if m is None:
            m = Z*(scipy.constants.m_p+scipy.constants.m_n)

        self.ions.append({'name': name, 'Z': int(Z), 'm': m})


    def getMass(self, name):
        """
        Returns the mass of the named ion species.
        """
        return self[name]['m']


    def __getitem__(self, name):
        for i in self.ions:
            if i['name'] == name:
                return i

        raise KeyError(f"No ion named '{name}' found among ions.")


    def getIndex(self, name):
        """
        Returns the index for the named ion.
        """
        for i in range(len(self.ions)):
            if self.ions[i]['name'] == name:
                return i

        raise Exception(f"No ion with name '{name}' found.")


    def __iter__(self):
        for ion in self.ions:
            yield ion


