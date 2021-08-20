# Load data from STREAM and initialize PYDYON with it

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
from . Equations import *
from . import ConfinementTime, IonHandler, PlasmaVolume, UnknownQuantityHandler


# Mapping of PYDYON terms to STREAM terms
# 
# PARAMETERS
#   pydyon:   Function constructing the PYDYON object to test.
#   stream:   Function extracting corresponding data to compare to from a STREAMOutput object.
# (optional)
#   eval:     Function to call to evaluate the PYDYON object. If not defined, called as 'obj(t,x)'.
TERMS = {
    'Radiated power': { 'pydyon': RadiatedPowerTerm, 'stream': lambda so : so.other.fluid.Tcold_radiation[:,0] },
    'Ohmic power': { 'pydyon': OhmicPowerTerm, 'stream': lambda so : -so.other.fluid.Tcold_ohmic[:,0] },
    'e-i equilibration': { 'pydyon': EquilibrationPowerTerm, 'stream': lambda so : so.other.fluid.Tcold_ion_coll[:,0] },
    'e heat convection': { 'pydyon': ElectronConvectivePowerTerm, 'stream': lambda so : so.other.scalar.energyloss_T_cold[:,0] },
    'Confinement time': { 'pydyon': ConfinementTime, 'stream': lambda so : so.other.stream.tau_D[:,0] },
    r'dI\_p / dt': { 'pydyon': CircuitEquation, 'eval': lambda ce, t, x : ce.dIp_dt(t,x), 'stream': lambda so : np.diff(so.eqsys.I_p[:,0]) / np.diff(so.grid.t[:]) },
    r'dI\_w / dt': { 'pydyon': CircuitEquation, 'eval': lambda ce, t, x : ce.dIMK2_dt(t,x), 'stream': lambda so : np.diff(so.eqsys.I_wall[:,0]) / np.diff(so.grid.t[:]) }
}

# Mapping from STREAMSettings to PYDYON settings.
SETTINGS = {
    'a': lambda ss : ss.radialgrid.a,
    'R': lambda ss : ss.radialgrid.R0,
    'V_vessel': lambda ss : ss.radialgrid.vessel_volume,
    'Bphi': lambda ss : ss.radialgrid.B0[0],
    'Bv': lambda ss : 1e-3,
    # The following currently only work with the 'TYPE_CIRCUIT' model in STREAM
    'l_MK2': lambda ss : ss.radialgrid.b,
    'Vloop': lambda ss : ss.eqsys.E_field.circuit_Vloop if ss.eqsys.E_field.circuit_Vloop.size==1 else scipy.interpolate.interp1d(ss.eqsys.E_field.circuit_Vloop_t, ss.eqsys.E_field.circuit_Vloop),
    'Lp': lambda ss : ss.eqsys.E_field.circuit_Lp,
    'LMK2': lambda ss : ss.eqsys.E_field.circuit_Lwall,
    'M': lambda ss : ss.eqsys.E_field.circuit_M,
    'RMK2': lambda ss : ss.eqsys.E_field.circuit_Rwall if ss.eqsys.E_field.circuit_Rwall<1 else np.inf
}


def compareToSTREAM(ss, so, verbose=True):
    # Load settings from STREAMSettings object
    settings = {}
    for s in SETTINGS.keys():
        settings[s] = SETTINGS[s](ss)

    # Construct ion object
    ions = IonHandler()
    ions.addIon('D', Z=1)

    pv = PlasmaVolume(a=settings['a'], R=settings['R'], V_vessel=settings['V_vessel'], ions=ions)

    # Construct unknown quantity handler
    unknowns = UnknownQuantityHandler(ions, pv)

    # Evaluate terms and plot those which deviate
    TOL = 1e-1
    for name, term in TERMS.items():
        t, pydyon, stream = _evaluateTerm(so, term, unknowns, ions, settings)

        Delta = np.amax(np.abs((pydyon-stream)/stream))

        if Delta > TOL:
            print("- Term '{}' is INACCURATE (delta = {:.3f}%)".format(name, Delta*100))

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(t, stream, 'k', label='STREAM')
            ax.plot(t, pydyon, 'r--', label='PYDYON')
            ax.set_xlabel('Time (s)')
            ax.set_title(name)
            ax.set_xlim([0, t[-1]])
            ax.legend()
            fig.tight_layout()
        else:
            print(f"- Term '{name}' ACCURATE to within 10%")

    plt.show()


def _evaluateTerm(so, term, unknowns, ions, settings):
    """
    Evaluate the term 'term', which can be either a key into the 'TERMS' dict
    or a dict defined in the same manner as the others appearing in 'TERMS'.

    :return: tuple consisting of time vector, PYDYON solution, STREAM solution.
    """
    if type(term) == str:
        term = TERMS[term]

    # Constructor or regular function?
    if '__init__' in term['pydyon'].__dict__:
        s = _getSubSettings(settings, term['pydyon'].__init__.__code__.co_varnames)
    else:
        s = _getSubSettings(settings, term['pydyon'].__code__.co_varnames)

    s['quantities'] = unknowns
    s['ions'] = ions

    # Instantiate PYDYON term
    obj = term['pydyon'](**s)
    # Get STREAM data
    stream = term['stream'](so)

    t = so.grid.t[1:]
    y = np.zeros((t.size,))
    for i in range(t.size):
        x = fromSTREAM(so, unknowns, time=i)

        if 'eval' in term:
            y[i] = term['eval'](obj, t[i], x)
        else:
            y[i] = obj(t[i], x)

    return t, y, stream


def _getSubSettings(allSettings, subset):
    """
    Extracts the keys named in 'subset' from the dictionary
    of all settings 'allSettings'.
    """
    s = {}
    for a in subset:
        if a in ['self', 'quantities', 'ions']:
            continue

        s[a] = allSettings[a]

    return s


def fromSTREAM(so, uqh, time=0, ion='D'):
    """
    Initialize a PYDYON unknown quantity vector with data from a
    STREAMOutput object.

    :param so:   STREAMOutput object to take data from.
    :param uqh:  PYDYON.UnknownQuantityHandler object to insert data with.
    :param time: Time index to initialize data from.
    :param ion:  Name of main ion species (to use for Ti).
    """
    dct = {
        'Te': so.eqsys.T_cold[time,0],
        'Ti': so.eqsys.W_i.getTemperature(ion)[time,0],
        'Ip': so.eqsys.I_p[time,0],
        'IMK2': so.eqsys.I_wall[time,0],
        # TODO insert all other ion densities
        f'ni{ion}': so.eqsys.n_i[ion].data[time,:,0]
    }

    return uqh.setvector(dct)

