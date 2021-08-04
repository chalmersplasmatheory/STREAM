# Representation of a STREAM mean-free path

from DREAM.Output.IonSpeciesFluidQuantity import IonSpeciesFluidQuantity


class MeanFreePath(IonSpeciesFluidQuantity):
    

    def __init__(self, *args, **kwargs):
        """
        Constructor.
        """
        super().__init__(*args, **kwargs)


