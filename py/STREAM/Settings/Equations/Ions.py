# Wrapper for handling ion species

import numpy as np
import DREAM.Settings.Equations.Ions as DREAMIons
from . IonSpecies import IonSpecies, IONS_PRESCRIBED, ION_OPACITY_MODE_TRANSPARENT


# Model to use for ion heat
IONS_T_I_NEGLECT = DREAMIons.IONS_T_I_NEGLECT
IONS_T_I_INCLUDE = DREAMIons.IONS_T_I_INCLUDE


class Ions(DREAMIons.Ions):
    

    def __init__(self, settings, *args, **kwargs):
        """
        Constructor.
        """
        super().__init__(settings=settings, *args, **kwargs)


    def addIon(self, name, Z, iontype=IONS_PRESCRIBED, Z0=None, isotope=0, SPIMolarFraction=-1, opacity_mode=ION_OPACITY_MODE_TRANSPARENT, T=None, n=None, r=None, t=None, tritium=False):
        """
        Adds a new ion species to the plasma.

        :param str name:               Name by which the ion species will be referred to.
        :param int Z:                  Ion charge number.
        :param int isotope:            Ion mass number.
        :param int iontype:            Method to use for evolving ions in time.
        :param int Z0:                 Charge state to populate (used for populating exactly one charge state for the ion).
        :param n:                      Ion density (can be either a scalar, 1D array or 2D array, depending on the other input parameters)
        :param float SPIMolarFraction: Molar fraction of the SPI injection (if any). A negative value means that this species is not part of the SPI injection 
        :param numpy.ndarray r:        Radial grid on which the input density is defined.
        :param T:                      Ion initial temperature (can be scalar for uniform temperature, otherwise 1D array matching `r` in size)
        :param numpy.ndarray r:        Radial grid on which the input density and temperature is defined.
        :param numpy.ndarray t:        Time grid on which the input density is defined.
        :param bool tritium:           If ``True``, the ion species is treated as Tritium.
        """
        if (self.r is not None) and (r is not None) and (np.any(self.r != r)):
            raise EquationException("The radial grid must be the same for all ion species.")
        if (self.t is not None) and (t is not None) and (np.any(self.t != t)):
            raise EquationException("The time grid must be the same for all ion species.")

        if T is not None:
            self.typeTi = IONS_T_I_INCLUDE

        ion = IonSpecies(settings=self.settings, name=name, Z=Z, ttype=iontype, Z0=Z0, isotope=isotope, SPIMolarFraction=SPIMolarFraction, opacity_mode=opacity_mode, T=T, n=n, r=r, t=t, interpr=self.r, interpt=None, tritium=tritium)

        self.ions.append(ion)

        self.r = ion.getR()
        if ion.getTime() is not None:
            self.t = ion.getTime()


    def setJET_CWrecycling(self, deuterium='D', carbon='C', oxygen='O', tritium='T', iron = 'Fe'):
        """
        Set recycling coefficients to the values estimated for JET-CW
        (as given in [Kim et al 2012 Nucl. Fusion 52 103016]).
        """
        for ion in self.ions:
            if ion.name == deuterium:
                ion.setRecyclingCoefficient(carbon, 0.015) #Ändrat från 0.03!!!
                ion.setRecyclingCoefficient(iron, 0.015)
            elif ion.name == oxygen:
                ion.setRecyclingCoefficient(carbon, 1.0)
                ion.setRecyclingCoefficient(oxygen, 1.0)
            elif ion.name == tritium:
                ion.setRecyclingCoefficient(tritium, 1.0)


    def fromdict(self, data):
        """
        Load ion data from the given dictionary.
        """
        super().fromdict(data)

        # Load recycling coefficients
        if 'recycling' in data:
            rec = data['recycling'][:]
            for i in range(len(self.ions)):
                for j in range(len(self.ions)):
                    if rec[i,j] != 0.0:
                        self.ions[i].setRecyclingCoefficient(self.ions[j].name, rec[i,j])


    def todict(self):
        """
        Convert this object into a Python dictionary.
        """
        data = super().todict()

        # Construct recycling coefficient table
        nions = len(self.ions)
        # recycling coefficient table
        rec = np.zeros((nions, nions))

        for i in range(nions):
            ion = self.ions[i]

            for j in range(nions):
                rec[i,j] = ion.getRecyclingCoefficient(self.ions[j].name)

        data['recycling'] = rec

        return data


