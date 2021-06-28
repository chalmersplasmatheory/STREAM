
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
        Retrieves the recycling coefficient value for the
        named ion species.
        """
        if species in self.recycling:
            return self.recycling[species]
        else:
            return 0.0


    def setRecyclingCoefficient(self, species, value):
        """
        Sets the recycling coefficient :math:`Y_i^j`, where :math:`i` denotes
        this ion species and :math:`j` the species of the interacting ion.

        :param str species: Name of interacting ion species.
        :param float value: Recycling coefficient value.
        """
        if species in self.recycling:
            print("WARNING: Updating recycling coefficient for species '{}' and '{}'.".format(self.name, species))

        self.recycling[species] = value


