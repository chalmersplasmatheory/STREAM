from STREAM.Settings.Equations.Ions import Ions
from DREAM.Settings.EquationSystem import EquationSystem as DREAMEqSys


class EquationSystem(DREAMEqSys):
    

    def __init__(self, settings):
        """
        Constructor.

        settings: Parent settings object.
        """
        super().__init__(settings)

        self.addUnknown('n_i', Ions(settings=settings))

        # TODO Overwrite with STREAM-modified unknowns

