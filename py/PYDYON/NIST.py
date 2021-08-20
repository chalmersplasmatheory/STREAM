# NIST data handler


import h5py
import numpy as np
import pathlib
import scipy.interpolate


class NIST:
    

    def __init__(self):
        """
        Constructor.
        """
        self.elements = {}

        path = str(pathlib.Path(__file__).parent.resolve().absolute())

        with h5py.File(path + '/NIST_ioniz.h5', 'r') as f:
            for el in f:
                self.load_element(el, f[el])

                if el == 'H':
                    self.load_element('D', f[el])


    def load_element(self, element, f):
        """
        Load data for an element using the given HDF5 file handle.
        """
        self.elements[element] = {
            'Z': int(f['Z'][:][0]),
            'data': f['data'][:]
        }


    def ionization(self, name, Z0):
        """
        Returns the ionization energy for the named ion species
        in the given charge state.
        """
        return self.elements[name]['data'][Z0]


