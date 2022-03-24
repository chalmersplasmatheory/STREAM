import numpy as np

import DREAM.Settings.Equations.IonSpecies as DREAMIonSpecies


# Types in DREAM
IONS_PRESCRIBED = DREAMIonSpecies.IONS_PRESCRIBED 
IONS_EQUILIBRIUM = DREAMIonSpecies.IONS_EQUILIBRIUM 
IONS_DYNAMIC = DREAMIonSpecies.IONS_DYNAMIC 

# Types which are extensions implemented in this interface
# (which are special cases of the DREAM types above)
IONS_DYNAMIC_NEUTRAL = DREAMIonSpecies.IONS_DYNAMIC_NEUTRAL 
IONS_DYNAMIC_FULLY_IONIZED = DREAMIonSpecies.IONS_DYNAMIC_FULLY_IONIZED 
IONS_PRESCRIBED_NEUTRAL = DREAMIonSpecies.IONS_PRESCRIBED_NEUTRAL 
IONS_PRESCRIBED_FULLY_IONIZED = DREAMIonSpecies.IONS_PRESCRIBED_FULLY_IONIZED 
IONS_EQUILIBRIUM_NEUTRAL = DREAMIonSpecies.IONS_EQUILIBRIUM_NEUTRAL 
IONS_EQUILIBRIUM_FULLY_IONIZED = DREAMIonSpecies.IONS_EQUILIBRIUM_FULLY_IONIZED

ION_OPACITY_MODE_TRANSPARENT = DREAMIonSpecies.ION_OPACITY_MODE_TRANSPARENT
ION_OPACITY_MODE_GROUND_STATE_OPAQUE = DREAMIonSpecies.ION_OPACITY_MODE_GROUND_STATE_OPAQUE

class IonSpecies(DREAMIonSpecies.IonSpecies):
    

    def __init__(self, *args, **kwargs):
        """
        Constructor.
        """
        super().__init__(*args, **kwargs)

        self.recycling = {}


    def getRecyclingCoefficient(self, species):
        """
        Retrieves the recycling coefficient value-vector and time vector for
        the named ion species.
        """
        if species in self.recycling:
            return self.recycling[species][0], self.recycling[species][1]
        else:
            return np.zeros(1), np.zeros(1)


    def pressureToDensity(p, T=300/11604.51812):
        """
        Convert a given prefill pressure (in Torr) to an atom density. The
        prefill pressure is related to the gas temperature T and density n via

          p = n*T
        
        Since 1 Torr = (1 atm)/760 = 133.32 Pa, we have that the atomic density
        at prefill pressure p (given in Torr) is

                            133.32
          (n [m^-3]) =  --------------- * (p [Torr])
                         k_B * (T [J])

        where k_B is Boltzmann's constant.

        :param p: Prefill gas pressure (in Torr).
        :param T: Temperature (in eV). Default: room temperature (= 300 K).
        """
        eV2J = 11604.51812
        return 133.32 * p / (scipy.constants.k * eV2J * T)


    def setRecyclingCoefficient(self, species, values, tvalues=[0.0]):
        """
        Sets the recycling coefficient :math:`Y_i^j`, where :math:`i` denotes
        this ion species and :math:`j` the species of the interacting ion.

        :param str species: Name of interacting ion species.
        :param float value: Recycling coefficient value.
        """
        if species in self.recycling:
            print("WARNING: Updating recycling coefficient for species '{}' and '{}'.".format(self.name, species))

        if isinstance(values, float) or isinstance(values, int):
            values = [float(values)]
        if isinstance(tvalues, float) or isinstance(tvalues, int):
            tvalues = [float(tvalues)]

        self.recycling[species] = [np.array(values), np.array(tvalues)]


