from DREAM.Settings.Equations.EquationException import EquationException
from DREAM.Settings.Equations.PrescribedParameter import PrescribedParameter
from DREAM.Settings.Equations.PrescribedInitialParameter import PrescribedInitialParameter
from DREAM.Settings.Equations.PrescribedScalarParameter import PrescribedScalarParameter
from DREAM.Settings.Equations.UnknownQuantity import UnknownQuantity

DEFAULT_TYPE = 1 #Bohm confinement time
DEFAULT_MIXED = 0 #Mixed confinement time is disabled
DEFAULT_SMOOTHLESS = 1 #Inverse sum


class ConfinementTime(UnknownQuantity):
    def __init__(self, settings, ttype=DEFAULT_TYPE, mmixed=DEFAULT_MIXED, ssmoothless = DEFAULT_SMOOTHLESS):
        super().__init__(settings=settings)
        self.type = ttype
        self.mixed = mmixed
        self.smoothless = ssmoothless

    def setType(self, ttype):
        self.type = ttype

    def setMixed(self, mmixed):
        self.mixed = mmixed
    
    def setSmoothless(self, ssmoothless) :
        self.smoothless = ssmoothless

    def todict(self):
        return dict(tau_perp=self.type, mixed=self.mixed, smoothless=self.smoothless)

    def fromdict(self, data):
        self.type = data['tau_perp']
        self.mixed = data['mixed']
        self.smoothless = data['smoothless']

    def verifySettings(self):
        if self.type < 1:
            raise EquationException("Invalid tau_perp type")
        else:
            pass
        if self.smoothless < 1:
            raise EquationException("Invalid smoothless type")
        else:
            pass

