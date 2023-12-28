from DREAM.Settings.Equations.EquationException import EquationException
from DREAM.Settings.Equations.PrescribedParameter import PrescribedParameter
from DREAM.Settings.Equations.PrescribedInitialParameter import PrescribedInitialParameter
from DREAM.Settings.Equations.PrescribedScalarParameter import PrescribedScalarParameter
from DREAM.Settings.Equations.UnknownQuantity import UnknownQuantity

DEFAULT_TYPE = 1


class ConfinementTime(UnknownQuantity):
    def __init__(self, settings, ttype=DEFAULT_TYPE):
        super().__init__(settings=settings)
        self.type = ttype

    def setType(self, ttype):
        self.type = ttype

    def todict(self):
        return dict(tau_perp=self.type)

    def fromdict(self, data):
        self.type = data['tau_perp']

    def verifySettings(self):
        if self.type < 1:
            raise EquationException("Invalid tau_perp type")
        else:
            pass
