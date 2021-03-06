
from DREAM.Settings.EquationSystem import EquationSystem as DREAMEqSys
from STREAM.Settings.Equations.ElectricField import ElectricField
from STREAM.Settings.Equations.Ions import Ions


class EquationSystem(DREAMEqSys):
    

    def __init__(self, settings):
        """
        Constructor.

        settings: Parent settings object.
        """
        super().__init__(settings)

        # Overwrite with STREAM-modified unknowns
        self.addUnknown('n_i', Ions(settings=settings))
        self.addUnknown('E_field', ElectricField(settings=settings))


