# Help class for accessing unknowns in the equation system

import numpy as np


class UnknownQuantityHandler:

    
    UNKNOWNS = [
        'We', 'Wi', 'Ip', 'IMK2'
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


    def __len__(self):
        """
        Returns the number of unknowns in the equation system.
        """
        return len(self.map)


    def __getitem__(self, name):
        """
        Returns data for the quantity with the given name.
        """
        if name.startswith('ni'):
            ion = self.ions[name[2:]]

            idx = self.map[f'ni{name[2:]}_0']
            return self.x[idx:(idx+ion['Z']+1)]
        elif name == 'ne':
            return self.getElectronDensity()
        elif name == 'Te':
            idx = self.map['We']
            We  = self.x[idx]
            ne  = self.getElectronDensity()
            return 2*We / (3*ne)
        elif name == 'Ti':
            idx = self.map['Wi']
            Wi  = self.x[idx]
            ni  = self.getTotalIonDensity()
            return 2*Wi / (3*ni)
        else:
            if name not in self.map:
                raise KeyError(f"No unknown quantity named '{name}' exists in equation system.")

            idx = self.map[name]

            return self.x[idx]


    def getElectronDensity(self):
        """
        Evaluates and returns the electron density.
        """
        ne = 0
        for ion in self.ions:
            A = ion['name']
            Z = ion['Z']
            ni = self.getIonData(A)
            
            for Z0 in range(1,Z0+1):
                ne += Z0 * ni[Z0]

        return ne


    def getIonData(self, ion):
        """
        Returns data for the named ion.
        """
        if type(ion) == dict:
            return self['ni{}'.format(ion['name'])]
        else:
            return self[f'ni{ion}']


    def getTotalIonDensity(self):
        """
        Calculates density of all ions (neutrals not included).
        """
        n = 0
        for ion in self.ions:
            ni = self.getIonData(ion)
            Z  = ion['Z']

            for Z0 in range(1, Z+1):
                n += ni[Z0]

        return n


    def getLambda(self, ion):
        """
        Returns the mean-free path for the named ion species.
        """
        ne = self['ne']
        Te = self['Te']
        Ti = self['Ti']

        return self.plasmavolume.getLambda(ion, ne=ne, Te=Te, Ti=Ti)


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


    def getZeff(self):
        """
        Calculates and returns the plasma effective charge.
        """
        Zup, Zdn = 0, 0
        for ion in self.ions:
            A  = ion['name']
            Z  = ion['Z']
            ni = self[f'ni{A}']

            for Z0 in range(1,Z+1):
                Zup += ni[Z0] * Z0**2
                Zdn += ni[Z0] * Z0

        return Zup / Zdn


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


