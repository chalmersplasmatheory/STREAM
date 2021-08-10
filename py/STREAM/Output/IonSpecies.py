# Specialized IonSpecies for STREAM

from DREAM.Output.IonSpecies import IonSpecies as DREAMIonSpecies
from . IonState import IonState


class IonSpecies(DREAMIonSpecies):
    

    def __init__(self, name, Z, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super().__init__(name=name, Z=Z, data=data, grid=grid, output=output, attr=attr)


    def addChargeState(self, name, Z, Z0, data, attr=list()):
        """
        Adds a new IonState object to the list of ion charge states.
        """
        self.ionstates.append(IonState(name=name, Z=Z, Z0=Z0, data=data, grid=self.grid, output=self.output, attr=attr))


