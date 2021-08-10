# Specialized IonHandler for STREAM

from DREAM.Output.IonHandler import IonHandler as DREAMIonHandler
from . IonSpecies import IonSpecies


class IonHandler(DREAMIonHandler):
    

    def __init__(self, name, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super().__init__(name=name, data=data, grid=grid, output=output, attr=attr)


    def addIon(self, name, Z, data, attr=list()):
        """
        Adds a new ion to the list of ions.
        """
        self.ions.append(IonSpecies(name=name, Z=Z, data=data, grid=self.grid, output=self.output, attr=attr))


