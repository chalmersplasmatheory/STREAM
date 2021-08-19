# Main equation solver of PYDYON


from scipy.integrate import solve_ivp
from . IonHandler import IonHandler
from . PlasmaVolume import PlasmaVolume
from . UnknownQuantityHandler import UnknownQuantityHandler

from . Equations import *


class Simulation:
    

    unknowns = None
    terms = []
    
    # Default settings
    settings = {
        # Tokamak parameters
        'a': 0.5,           # Minor radius (m)
        'R': 3,             # Major radius (m)
        'V_vessel': 100,    # Vessel volume (m^3)
        'Bphi': 2.3,        # Toroidal magnetic field (T)
        'Bv': 1e-3,         # Vertical magnetic field (T)
        'l_MK2': 1,         # Distance between plasma centre and passive structure (m)
    }


    def __init__(self, **settings):
        """
        Constructor.

        :param settings: Dictionary containing simulation parameters.
        """
        self.ions = IonHandler()


    def addIon(self, name, Z, m=None):
        """
        Add ions to the simulation.

        :param name: Name of ion.
        :param Z:    Atomic charge number.
        :param m:    Atomic mass (if ``None``, calculated as Z*(m_p + m_n)).
        """
        self.ions.addIon(name=name, Z=Z, m=m)


    def _constructEquationSystem(self):
        """
        Constructs the equation system.
        """
        def i(name): return self.unknowns.map(name)
        eqsys = [None] * len(self.unknowns)

        s = self.settings

        #######
        # Power balance equations
        tausettings = {'Bphi': s['Bphi'], 'Bv': s['Bv'], 'l_MK2': s['l_MK2']}
        Poh    = OhmicPowerTerm(self.unknowns, self.ions)
        Prad   = RadiatedPowerTerm(self.unknowns, self.ions)
        Pequi  = EquilibrationPowerTerm(self.unknowns, self.ions)
        Pconve = ElectronConvectivePowerTerm(self.unknowns, self.ions, **tausettings)
        Pconvi = IonConvectivePowerTerm(self.unknowns, self.ions, **tausettings)
        Pcx    = ChargeExchangePowerTerm(self.unknowns, self.ions)

        # Ensure that terms are not deleted by the garbage collector
        self.terms.extend((Poh, Prad, Pequi, Pconve, Pconvi, Pcx))

        eqsys[i('We')] = lambda t, x : Poh(x) - Prad(x) - Pequi(x) - Pconve(x)
        eqsys[i('Wi')] = lambda t, x : Pequi(x) - Pcx(x) - Pconvi(x)

        ######
        # Circuit equation
        circuit = CircuitEquation(


    def initialize(self, **values):
        """
        Set the initial values of all quantities.
        """
        if len(self.ions) == 0:
            raise Exception("Ions must be added before calling 'initialize()'.")

        s = self.updateSettings(self.settings, values)

        pv = PlasmaVolume(a=s['a'], R=s['R'], V_vessel=s['V_vessel'], ions=self.ions)
        self.unknowns = UnknownQuantityHandler(self.ions, pv)
        self.unknowns.setvector(values)


    def solve(uqh):
        """
        Solve the system of equations using the given settings.

        :param uqh: UnknownQuantityHandler of 
        """
        dydt = self._constructEquationSystem()
        pass


    def updateSettings(self, set1, set2):
        """
        Update the settings defined in the dict 'set1' with values given
        in the dict 'set2'. Entries which are only defined in 'set1' are
        kept as is, while entries which only exist in 'set2' raise an
        exception.
        """
        for k, v in set2.items():
            if k not in set1:
                raise KeyError(f"No setting '{k}' provided by PYDYON.")
            else:
                set1[k] = v

        return set1

