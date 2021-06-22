
from DREAM.Settings.Equations.IonSpecies import IonSpecies as DREAMIonSpecies


class IonSpecies(DREAMIonSpecies):
    

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


