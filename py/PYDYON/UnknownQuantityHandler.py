# Help class for accessing unknowns in the equation system

import numpy as np
import scipy.constants


class UnknownQuantityHandler:

    
    UNKNOWNS = [
        'We', 'Wi', 'Ip', 'IMK2'
    ]
    

    def __init__(self, ions, plasmavolume):
        """
        Constructor.
        """
        self.t = None
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
        e = scipy.constants.e

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
            return 2*We / (3*e*ne)
        elif name == 'Ti':
            idx = self.map['Wi']
            Wi  = self.x[idx]
            ni  = self.getTotalIonDensity()
            return 2*Wi / (3*e*ni)
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
            
            for Z0 in range(1,Z+1):
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


    def getNameByIndex(self, index):
        """
        Returns the name of the unknown quantity with
        the given index.
        """
        for u, i in self.map.items():
            if index == i:
                return u

        raise KeyError(f"No unknown with index {index}.")


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


    def getV_n(self, ion): return self.getNeutralVolume(ion)


    def getV_n_tot(self, ion): return self.getTotalNeutralVolume(ion)


    def getV_p(self): return self.getPlasmaVolume()


    def getNeutralVolume(self, ion):
        """
        Returns the neutral volume for the named ion species.
        """
        ne = self['ne']
        Te = self['Te']
        Ti = self['Ti']

        return self.plasmavolume.getV_n(self.t, ion, ne=ne, Te=Te, Ti=Ti)
    

    def getTotalNeutralVolume(self, ion):
        """
        Returns the total neutral volume for the named ion species.
        """
        ne = self['ne']
        Te = self['Te']
        Ti = self['Ti']

        return self.plasmavolume.getV_n_tot(self.t, ion, ne=ne, Te=Te, Ti=Ti)


    def getPlasmaVolume(self):
        """
        Returns the plasma volume V_p.
        """
        return self.plasmavolume.getV_p(self.t)


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


    def update(self, t, x):
        """
        Set new data for the unknown quantities.
        """
        self.t = t
        self.x = x


    def getdict(self, x=None):
        """
        Convert a solution vector to a dictionary.
        """
        e = scipy.constants.e

        if x is None:
            x = self.x

        dct = {}
        for k, idx in self.map.items():
            if k.startswith('ni'):
                nm, _, Z0 = k.partition('_')
                Z0 = int(Z0)
                A = nm[2:]
                if nm not in dct:
                    ion = self.ions[A]
                    dct[nm] = np.zeros((ion['Z']+1, len(x[idx])))
                dct[nm][Z0,:] = x[idx]
            else:
                dct[k] = x[idx]

        # Additional quantities
        ne = np.zeros((dct['We'].size, ))
        ni = np.zeros(ne.shape)
        for ion in self.ions:
            A = ion['name']
            Z = ion['Z']
            dens = dct[f'ni{A}']
            for Z0 in range(1,Z+1):
                ne += Z0*dens[Z0,:]
                ni += dens[Z0,:]

        dct['ne'] = ne
        dct['Te'] = (2.0/3.0) * dct['We'] / (e*ne)
        dct['Ti'] = (2.0/3.0) * dct['Wi'] / (e*ni)

        return dct


    def setvector(self, dct=None, t=0, **kwargs):
        """
        Converts a dictionary with names for unknowns to a numpy
        array.
        """
        e = scipy.constants.e

        if dct is None:
            dct = kwargs

        x = np.zeros((len(self.map),))
        Te, Ti = None, None
        
        for k, v in dct.items():
            # We set ion values in bulk
            if k.startswith('ni'):
                ion = self.ions[k[2:]]
                idx = self.map[f'ni{k[2:]}_0']
                x[idx:(idx+ion['Z']+1)] = dct[k]
            elif k == 'Te':
                # Save for later...
                Te = dct[k]
            elif k == 'Ti':
                # Save for later...
                Ti = dct[k]
            else:
                x[self.map[k]] = v

        self.update(t=t, x=x)

        if Te is not None:
            x[self.map['We']] = 3/2 * e*Te * self['ne']
        if Ti is not None:
            x[self.map['Wi']] = 3/2 * e*Ti * self.getTotalIonDensity()

        self.update(t=t, x=x)

        return x


