from DREAM.Settings.Equations.EquationException import EquationException
from DREAM.Settings.Equations.PrescribedParameter import PrescribedParameter
from DREAM.Settings.Equations.PrescribedInitialParameter import PrescribedInitialParameter
from DREAM.Settings.Equations.PrescribedScalarParameter import PrescribedScalarParameter
from DREAM.Settings.Equations.UnknownQuantity import UnknownQuantity

DEFAULT_TYPE = 1 #Bohm confinement time
DEFAULT_MIXED = 0 #Mixed confinement time is disabled


class ConfinementTime(UnknownQuantity):
    def __init__(self, settings, ttype=DEFAULT_TYPE, mmixed=DEFAULT_MIXED):
        super().__init__(settings=settings)
        self.type = ttype
        self.mixed = mmixed

    def setType(self, ttype):
        self.type = ttype

    def setMixed(self, mmixed):
        self.mixed = mmixed

    def todict(self):
        return dict(tau_perp=self.type, mixed=self.mixed)

    def fromdict(self, data):
        self.type = data['tau_perp']
        self.mixed = data['mixed']

    def verifySettings(self):
        if self.type < 1:
            raise EquationException("Invalid tau_perp type")
        else:
            pass

