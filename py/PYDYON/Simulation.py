# Main equation solver of PYDYON


from scipy.integrate import solve_ivp
from . IonHandler import IonHandler
from . PlasmaVolume import PlasmaVolume
from . SimulationResult import SimulationResult
from . UnknownQuantityHandler import UnknownQuantityHandler

from . Equations import *


class Simulation:
    

    unknowns = None
    _terms = []
    
    # Default settings
    settings = {
        # Tokamak parameters
        'a': 0.5,           # Minor radius (m)
        'R': 3,             # Major radius (m)
        'V_vessel': 100,    # Vessel volume (m^3)
        'Bphi': 2.3,        # Toroidal magnetic field (T)
        'Bv': 1e-3,         # Vertical magnetic field (T)
        'l_MK2': 1,         # Distance between plasma centre and passive structure (m)

        # Electrical parameters
        'Vloop': 20,        # External loop voltage (V)
        'Lp': 5.4e-6,       # Plasma inductance (H)
        'LMK2': 9.1e-6,     # Inductance of passive structure (H)
        'M': 2.49e-6,       # Mutual inductance of plasma and passive structure (H)
        'RMK2': 7.5e-4      # Resistance of passive structure (Ohm)
    }


    def __init__(self, **settings):
        """
        Constructor.

        :param settings: Dictionary containing simulation parameters.
        """
        self.ions = IonHandler()

        s = self.updateSettings(self.settings, settings)


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
        def i(name): return self.unknowns.map[name]
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
        self._terms.extend((Poh, Prad, Pequi, Pconve, Pconvi, Pcx))

        eqsys[i('We')] = lambda t, x : Poh(t, x) - Prad(t, x) - Pequi(t, x) - Pconve(t, x)
        eqsys[i('Wi')] = lambda t, x : Pequi(t, x) - Pcx(t, x) - Pconvi(t, x)

        #######
        # Circuit equation
        circuitsettings = {
            'Vloop': s['Vloop'],
            'Lp':    s['Lp'],
            'LMK2':  s['LMK2'],
            'M':     s['M'],
            'RMK2':  s['RMK2']
        }
        circuit = CircuitEquation(self.unknowns, self.ions, **circuitsettings)
        self._terms.append(circuit)

        eqsys[i('Ip')]   = lambda t, x : circuit.dIp_dt(t, x)
        eqsys[i('IMK2')] = lambda t, x : circuit.dIMK2_dt(t, x)

        #######
        # Ion equations
        Dizcx0  = DeuteriumAtomBalance(self.unknowns, self.ions)
        Dizcx1  = DeuteriumIonBalance(self.unknowns, self.ions)
        Din     = DeuteriumInflux(self.unknowns, self.ions, simple=True, **tausettings)
        itransp = IonTransport(self.unknowns, self.ions, **tausettings)
        Iizcx   = IonParticleBalance(self.unknowns, self.ions)
        Iin     = IonInflux(self.unknowns, self.ions, **tausettings)

        self._terms.extend((Dizcx0, Dizcx1, Din, itransp, Iizcx, Iin))
        for ion in self.ions:
            A = ion['name']
            Z = ion['Z']

            if A == 'D':
                eqsys[i('niD_0')] = lambda t, x : Dizcx0(t, x) + Din(t, x)
                eqsys[i('niD_1')] = lambda t, x : Dizcx1(t, x) - itransp(t, x, 'D', Z0=1)
            else:
                # Add neutral equations
                for Z0 in range(1, Z+1):
                    if Z0 == 0:
                        eqsys[i(f'ni{A}_{Z0}')] = lambda t, x : Iizcx(t, x, A, Z0) + Iin(t, x, A)
                    else:
                        eqsys[i(f'ni{A}_{Z0}')] = lambda t, x : Iizcx(t, x, A, Z0) - Itransp(t, x, A, Z0=Z0)

        atol = [0.0] * len(eqsys)

        atol[i('Ip')] = 1e-1
        atol[i('IMK2')] = 1

        return eqsys, atol


    def initialize(self, **values):
        """
        Set the initial values of all quantities.
        """
        if len(self.ions) == 0:
            raise Exception("Ions must be added before calling 'initialize()'.")

        s = self.settings

        pv = PlasmaVolume(a=s['a'], R=s['R'], V_vessel=s['V_vessel'], ions=self.ions)
        self.unknowns = UnknownQuantityHandler(self.ions, pv)
        self.unknowns.setvector(values)


    def solve(self, tMax):
        """
        Solve the system of equations using the given settings.

        :param tMax: Duration of simulation.
        """
        equations, atol = self._constructEquationSystem()

        def dydt(t, x):
            print(f't = {t} s')
            self.unknowns.update(x)
            return [f(t,x) for f in equations]

        sol = solve_ivp(dydt, t_span=(0, tMax), y0=self.unknowns.x, method='Radau')
        
        return SimulationResult(sol.t, self.unknowns.getdict(x=sol.y), simulation=self)


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

