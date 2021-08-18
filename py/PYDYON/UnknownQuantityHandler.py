# Help class for accessing unknowns in the equation system

import numpy as np


class UnknownQuantityHandler:

    
    UNKNOWNS = [
        'ne', 'Te', 'Ti'
    ]
    

    def __init__(self, ions, plasmavolume):
        """
        Constructor.
        """
        self.x = None
        self.plasmavolume = plasmavolume
        self.ions = ions

        self.map = {}
        idx = 0
        for u in self.UNKNOWNS:
            self.map[u] = idx
            idx += 1

        # Add ions
        for ion in self.ions:
            nm, Z = ion['name'], ion['Z']
            
            for i in range(Z+1):
                self.map[f'ni{nm}_{i}'] = idx
                idx += 1


    def __getitem__(self, name):
        """
        Returns data for the quantity with the given name.
        """
        if name.startswith('ni'):
            ion = self.ions[name[2:]]

            idx = self.map[f'ni{name[2:]}_0']
            return self.x[idx:(idx+ion['Z']+1)]
        else:
            if name not in self.map:
                raise KeyError(f"No unknown quantity named '{name}' exists in equation system.")

            idx = self.map[name]

            return self.x[idx]


    def getIonData(self, ion):
        """
        Returns data for the named ion.
        """
        if type(ion) == dict:
            return self['ni{}'.format(ion['name'])]
        else:
            return self[f'ni{ion}']


    def getNeutralVolume(self, ion):
        """
        Returns the neutral volume for the named ion species.
        """
        ne = self['ne']
        Te = self['Te']
        Ti = self['Ti']

        return self.plasmavolume.getV_n(ion, ne=ne, Te=Te, Ti=Ti)


    def getPlasmaVolume(self):
        """
        Returns the plasma volume V_p.
        """
        return self.plasmavolume.getV_p()


    def update(self, x):
        """
        Set new data for the unknown quantities.
        """
        self.x = x


    def setvector(self, dct):
        """
        Converts a dictionary with names for unknowns to a numpy
        array.
        """
        x = np.zeros((len(self.map),))
        
        for k, v in dct.items():
            # We set ion values in bulk
            if k.startswith('ni'):
                ion = self.ions[k[2:]]
                idx = self.map[f'ni{k[2:]}_0']
                x[idx:(idx+ion['Z']+1)] = dct[k]
            else:
                x[self.map[k]] = v

        self.update(x)
        return x


