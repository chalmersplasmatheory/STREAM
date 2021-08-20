# ADAS data handler


import h5py
import numpy as np
import pathlib
import scipy.interpolate


class ADAS:
    

    def __init__(self):
        """
        Constructor.
        """
        self.ACD = ADASRate('ACD', shiftup=True)    # Recombination
        self.CCD = ADASRate('CCD', shiftup=True)    # Charge-exchange
        self.PLT = ADASRate('PLT')                  # Line radiation
        self.PRB = ADASRate('PRB', shiftup=True)    # Recombination + bremsstrahling radiation
        self.SCD = ADASRate('SCD')                  # Ionization

        path = str(pathlib.Path(__file__).parent.resolve().absolute())

        with h5py.File(path + '/ADAS.h5', 'r') as f:
            for el in f:
                self.load_element(el, f[el])


    def load_element(self, element, f):
        """
        Load data for an element using the given HDF5 file handle.
        """
        self.ACD.addElement(element, f)
        self.CCD.addElement(element, f)
        self.PLT.addElement(element, f)
        self.PRB.addElement(element, f)
        self.SCD.addElement(element, f)


class ADASRate:
    

    def __init__(self, name, shiftup=False):
        """
        ADAS rate object.
        """
        self.elements = {}
        self.name = name
        self.shiftup = shiftup

    
    def addElement(self, name, f):
        """
        Extract data for this ADAS rate for the named element
        from the given HDF5 handle.
        """
        x = f[self.name.lower()]['data'][:]
        n = f[self.name.lower()]['n'][:]
        T = f[self.name.lower()]['T'][:]

        A = int(f['A'][:][0])
        Z = int(f['Z'][:][0])

        interp = []
        if self.shiftup:
            interp.append(None)

        for i in range(x.shape[0]):
            interp.append(scipy.interpolate.interp2d(n, T, x[i,:], kind='linear', bounds_error=False))

        # Make sure 'interp' has Z+1 elements
        if len(interp) < Z+1:
            interp.append(None)

        self.elements[name] = {
            'A': A,
            'Z': Z,
            'rates': interp
        }


    def __call__(self, *args, **kwargs):
        return self.eval(*args, **kwargs)


    def eval(self, name, Z0, n, T):
        """
        Evaluate this rate for the named element and charge state.
        """
        r = self.elements[name]['rates'][Z0]

        if r is None:
            raise ValueError(f"Cannot evaluate ADAS rate {self.name} for Z0 = {Z0}.")

        ln, lT = np.log10(n), np.log10(T)
        exp = r(ln, lT)

        if exp.size == 1:
            return 10**exp[0]
        else:
            return 10**exp


