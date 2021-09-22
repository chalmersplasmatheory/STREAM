# Load data from STREAM and initialize PYDYON with it

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
from . Equations import *
from . import ConfinementTime, IonHandler, PlasmaVolume, UnknownQuantityHandler, Simulation


# Mapping of PYDYON terms to STREAM terms
# 
# PARAMETERS
#   pydyon:   Function constructing the PYDYON object to test.
#   stream:   Function extracting corresponding data to compare to from a STREAMOutput object.
# (optional)
#   atol:     Absolute tolerance required for being considered accurate (default: 0).
#   eval:     Function to call to evaluate the PYDYON object. If not defined, called as 'obj(t,x)'.
#   
TERMS = {
    # Volume terms
    #'Plasma volume': {'pydyon': PlasmaVolume, 'eval': lambda pv, t, x, _ : pv.getV_p(t),'stream': lambda so: so.other.stream.V_p[:].flatten()},
    #'D neutral volume': {'pydyon': PlasmaVolume, 'eval': lambda pv, t, _, uqh : pv.getV_n(t, 'D', uqh['ne'], uqh['Te'], uqh['Ti']),'stream': lambda so: so.other.stream.V_n['D'][:].flatten()},
    #'C neutral volume': {'pydyon': PlasmaVolume, 'eval': lambda pv, t, _, uqh : pv.getV_n(t, 'C', uqh['ne'], uqh['Te'], uqh['Ti']),'stream': lambda so: so.other.stream.V_n['C'][:].flatten()},
    #'O neutral volume': {'pydyon': PlasmaVolume, 'eval': lambda pv, t, _, uqh : pv.getV_n(t, 'O', uqh['ne'], uqh['Te'], uqh['Ti']),'stream': lambda so: so.other.stream.V_n['O'][:].flatten()},
    #'D total neutral volume': {'pydyon': PlasmaVolume, 'eval': lambda pv, t, _, uqh : pv.getV_n_tot(t, 'D', uqh['ne'], uqh['Te'], uqh['Ti']),'stream': lambda so: so.other.stream.V_n_tot['D'][:].flatten()},
    #'C total neutral volume': {'pydyon': PlasmaVolume, 'eval': lambda pv, t, _, uqh : pv.getV_n_tot(t, 'C', uqh['ne'], uqh['Te'], uqh['Ti']),'stream': lambda so: so.other.stream.V_n_tot['C'][:].flatten()},
    #'O total neutral volume': {'pydyon': PlasmaVolume, 'eval': lambda pv, t, _, uqh : pv.getV_n_tot(t, 'O', uqh['ne'], uqh['Te'], uqh['Ti']),'stream': lambda so: so.other.stream.V_n_tot['O'][:].flatten()},
    #'LambdaD': {'pydyon': PlasmaVolume, 'eval': lambda pv, t, _, uqh : pv.getLambda('D', uqh['ne'], uqh['Te'], uqh['Ti']),'stream': lambda so: so.eqsys.lambda_i['D'][1:,0]},
    #'LambdaC': {'pydyon': PlasmaVolume, 'eval': lambda pv, t, _, uqh: pv.getLambda('C', uqh['ne'], uqh['Te'], uqh['Ti']), 'stream': lambda so: so.eqsys.lambda_i['C'][1:,0]},
    #'LambdaO': {'pydyon': PlasmaVolume, 'eval': lambda pv, t, _, uqh: pv.getLambda('O', uqh['ne'], uqh['Te'], uqh['Ti']), 'stream': lambda so: so.eqsys.lambda_i['O'][1:,0]},

    #'Radiated power': { 'pydyon': RadiatedPowerTerm, 'stream': lambda so : so.other.fluid.Tcold_radiation[:,0] },
    #'D Ion temperature': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalTi(t, x), 'stream': lambda so : 2.0/3.0*so.eqsys.W_i['D'][:-1]/so.eqsys.N_i['D'][:-1]/1.60217662e-19 },
    #'C Ion temperature': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalTi(t, x), 'stream': lambda so : 2.0/3.0*so.eqsys.W_i['C'][:-1]/so.eqsys.N_i['C'][:-1]/1.60217662e-19 },
    #'O Ion temperature': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalTi(t, x), 'stream': lambda so : 2.0/3.0*so.eqsys.W_i['O'][:-1]/so.eqsys.N_i['O'][:-1]/1.60217662e-19 },

    #'D-0 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniD(t, x)[0], 'stream': lambda so : so.eqsys.n_i['D'][0][1:]},
    #'D-1 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniD(t, x)[1], 'stream': lambda so : so.eqsys.n_i['D'][1][1:]},
    #'C-0 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniC(t, x)[0], 'stream': lambda so : so.eqsys.n_i['C'][0][1:]},
    #'C-1 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniC(t, x)[1], 'stream': lambda so : so.eqsys.n_i['C'][1][1:]},
    #'C-2 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniC(t, x)[2], 'stream': lambda so : so.eqsys.n_i['C'][2][1:]},
    #'C-3 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniC(t, x)[3], 'stream': lambda so : so.eqsys.n_i['C'][3][1:]},
    #'C-4 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniC(t, x)[4], 'stream': lambda so : so.eqsys.n_i['C'][4][1:]},
    #'C-5 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniC(t, x)[5], 'stream': lambda so : so.eqsys.n_i['C'][5][1:]},
    #'C-6 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniC(t, x)[6], 'stream': lambda so : so.eqsys.n_i['C'][6][1:]},
    #'O-0 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniO(t, x)[0], 'stream': lambda so : so.eqsys.n_i['O'][0][1:]},
    #'O-1 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniO(t, x)[1], 'stream': lambda so : so.eqsys.n_i['O'][1][1:]},
    #'O-2 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniO(t, x)[2], 'stream': lambda so : so.eqsys.n_i['O'][2][1:]},
    #'O-3 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniO(t, x)[3], 'stream': lambda so : so.eqsys.n_i['O'][3][1:]},
    #'O-4 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniO(t, x)[4], 'stream': lambda so : so.eqsys.n_i['O'][4][1:]},
    #'O-5 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniO(t, x)[5], 'stream': lambda so : so.eqsys.n_i['O'][5][1:]},
    #'O-6 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniO(t, x)[6], 'stream': lambda so : so.eqsys.n_i['O'][6][1:]},
    #'O-7 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniO(t, x)[7], 'stream': lambda so : so.eqsys.n_i['O'][7][1:]},
    #'O-8 density': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.evalniO(t, x)[8], 'stream': lambda so : so.eqsys.n_i['O'][8][1:]},

    'Ohmic power': { 'pydyon': OhmicPowerTerm, 'stream': lambda so : -so.other.fluid.Tcold_ohmic[:,0] },
    'e-i equilibration': { 'pydyon': EquilibrationPowerTerm, 'stream': lambda so : so.other.fluid.Tcold_ion_coll[:,0] },
    'i-e equilibration': { 'pydyon': EquilibrationPowerTerm, 'stream': lambda so : so.other.stream.Wi_e_coll[:,0] },
    #'e heat convection': { 'pydyon': ElectronConvectivePowerTerm, 'stream': lambda so : so.other.scalar.energyloss_T_cold[:,0] },
    #'e heat convection': { 'pydyon': ElectronConvectivePowerTerm, 'stream': lambda so : so.other.stream.Tcold_transport[:,0] },
    'i heat convection': { 'pydyon': IonConvectivePowerTerm, 'stream': lambda so : -np.sum(so.other.stream.Wi_iontransport[:,:,0], axis=1) },
    #'Charge-exchange heat loss': { 'pydyon': ChargeExchangePowerTerm, 'stream': lambda so : -so.other.stream.Wi_chargeexchange[:,0] },
    #'i particle transport': { 'pydyon': IonTransport, 'eval': lambda ce, t, x, _: ce(t, x, 'D', Z0=1), 'stream': lambda so : -so.other.stream.ni_iontransport[:,1,0] },
    #'Confinement time': { 'pydyon': ConfinementTime, 'stream': lambda so : so.other.stream.tau_D[:,0] },
    #r'dI\_p / dt': { 'pydyon': CircuitEquation, 'eval': lambda ce, t, x, _: ce.dIp_dt(t,x), 'stream': lambda so : np.diff(so.eqsys.I_p[:,0]) / np.diff(so.grid.t[:]) , 'atol': 1 },
    #r'dI\_w / dt': { 'pydyon': CircuitEquation, 'eval': lambda ce, t, x, _ : ce.dIMK2_dt(t,x), 'stream': lambda so : np.diff(so.eqsys.I_wall[:,0]) / np.diff(so.grid.t[:]), 'atol': 1 },

    # Specialized deuterium ionization
    #'D-0 positive ionization': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.eval(t, x, True)[1], 'stream': lambda so : so.other.stream.ionrateequation_posIonization[:,0,0] },
    #'D-0 negative ionization': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.eval(t, x, True)[2], 'stream': lambda so : so.other.stream.ionrateequation_negIonization[:,0,0] },
    #'D-0 positive recombination': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.eval(t, x, True)[3], 'stream': lambda so : so.other.stream.ionrateequation_posRecombination[:,0,0] },
    #'D-0 negative recombination': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.eval(t, x, True)[4], 'stream': lambda so : so.other.stream.ionrateequation_negRecombination[:,0,0] },
    #'D-0 positive C-X': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.eval(t, x, True)[5], 'stream': lambda so : so.other.stream.ionrateequation_posChargeExchange[:,0,0] },
    #'D-0 negative C-X': { 'pydyon': DeuteriumAtomBalance, 'eval': lambda dab, t, x, _ : dab.eval(t, x, True)[6], 'stream': lambda so : so.other.stream.ionrateequation_negChargeExchange[:,0,0] },
    #'D-1 positive ionization': { 'pydyon': DeuteriumIonBalance, 'eval': lambda dab, t, x, _ : dab.eval(t, x, True)[1], 'stream': lambda so : so.other.stream.ionrateequation_posIonization[:,1,0] },
    #'D-1 negative ionization': { 'pydyon': DeuteriumIonBalance, 'eval': lambda dab, t, x, _ : dab.eval(t, x, True)[2], 'stream': lambda so : so.other.stream.ionrateequation_negIonization[:,1,0] },
    #'D-1 positive recombination': { 'pydyon': DeuteriumIonBalance, 'eval': lambda dab, t, x, _ : dab.eval(t, x, True)[3], 'stream': lambda so : so.other.stream.ionrateequation_posRecombination[:,1,0] },
    #'D-1 negative recombination': { 'pydyon': DeuteriumIonBalance, 'eval': lambda dab, t, x, _ : dab.eval(t, x, True)[4], 'stream': lambda so : so.other.stream.ionrateequation_negRecombination[:,1,0] },
    #'D-1 positive C-X': { 'pydyon': DeuteriumIonBalance, 'eval': lambda dab, t, x, _ : dab.eval(t, x, True)[5], 'stream': lambda so : so.other.stream.ionrateequation_posChargeExchange[:,1,0] },
    #'D-1 negative C-X': { 'pydyon': DeuteriumIonBalance, 'eval': lambda dab, t, x, _ : dab.eval(t, x, True)[6], 'stream': lambda so : so.other.stream.ionrateequation_negChargeExchange[:,1,0] },
    #'D-0 influx': { 'pydyon': DeuteriumInflux, 'stream': lambda so : so.other.stream.neutralinflux['D'][:] / so.other.stream.V_n_tot['D'][:] }
}

# Mapping from STREAMSettings to PYDYON settings.
SETTINGS = {
    'a': lambda ss : ss.radialgrid.a,
    'ta': lambda ss : ss.radialgrid.ta,
    'R': lambda ss : ss.radialgrid.R0,
    'V_vessel': lambda ss : ss.radialgrid.vessel_volume,
    'kappa': lambda ss : ss.radialgrid.kappa,
    'Bphi': lambda ss : ss.radialgrid.B0[0],
    'Bv': lambda ss : 1e-3,
    # The following currently only work with the 'TYPE_CIRCUIT' model in STREAM
    'l_MK2': lambda ss : ss.radialgrid.b,
    'Vloop': lambda ss : ss.eqsys.E_field.circuit_Vloop if ss.eqsys.E_field.circuit_Vloop.size==1 else scipy.interpolate.interp1d(ss.eqsys.E_field.circuit_Vloop_t, ss.eqsys.E_field.circuit_Vloop),
    'Lp': lambda ss : ss.eqsys.E_field.circuit_Lp,
    'LMK2': lambda ss : ss.eqsys.E_field.circuit_Lwall,
    'M': lambda ss : ss.eqsys.E_field.circuit_M,
    'RMK2': lambda ss : ss.eqsys.E_field.circuit_Rwall if ss.eqsys.E_field.circuit_Rwall<1 else np.inf,
    # "simple" neutral deuterium influx
    'simple': lambda ss : ss.radialgrid.c2==1
}


def compareToSTREAM(ss, so, verbose=True):
    """
    Compare individual terms in the equations solved by PYDON
    to their STREAM equivalents.
    """
    settings = loadSTREAMSettings(ss)

    # Construct ion object
    ions = IonHandler()
    for ion in so.eqsys.n_i.ions:
        ions.addIon(ion.name,ion.Z)
    #ions.addIon('D', Z=1)
    #ions.addIon('C', Z=6)
    #ions.addIon('O', Z=8)

    pv = PlasmaVolume(a=settings['a'], R=settings['R'], V_vessel=settings['V_vessel'], ions=ions, ta=settings['ta'])

    # Construct unknown quantity handler
    unknowns = UnknownQuantityHandler(ions, pv)

    # Evaluate terms and plot those which deviate
    RTOL = 1e-1
    for name, term in TERMS.items():
        t, pydyon, stream = _evaluateTerm(so, term, unknowns, ions, settings)
        ATOL = 0 if 'atol' not in term else term['atol']

        Delta = np.abs(pydyon-stream)

        if np.any(Delta > (RTOL*np.abs(stream) + ATOL)):
            if verbose:
                print("- Term '{}' is INACCURATE (delta = {:.3f}%)".format(name, np.amax(np.abs(Delta/stream))*100))

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(t, stream, 'k', label='STREAM')
            ax.plot(t, pydyon, 'r--', label='PYDYON')
            #ax.plot(t, stream/pydyon, 'k')
            ax.set_xlabel('Time (s)')
            ax.set_title(name)
            ax.set_xlim([0, t[-1]])
            ax.legend()
            fig.tight_layout()
        else:
            if verbose:
                print(f"- Term '{name}' ACCURATE to within 10%")

    plt.show()


def compareToSTREAMdt(ss, so, verbose=True):
    """
    Compare the time derivatives of unknown quantities solved for
    in DYON to their STREAM equivalents.
    """
    settings = loadSTREAMSettings(ss)
    if 'simple' in settings: del settings['simple']

    sim = Simulation(**settings)
    sim.addIon('D', 1)

    nD0 = np.array([so.eqsys.n_i['D'][0][0,0], so.eqsys.n_i['D'][1][0,0]])
    sim.initialize(We=so.eqsys.W_cold[0,0], Wi=so.eqsys.W_i['D'][0,0],
                   Ip=so.eqsys.I_p[0,0], IMK2=so.eqsys.I_wall[0,0],
                   niD=nD0)
    uqh = sim.unknowns

    equations, _ = sim._constructEquationSystem()

    dWedt = []
    dWidt = []
    dIpdt = []
    dIMK2dt = []
    dnD0dt = []
    dnD1dt = []

    t = so.grid.t
    sdWedt = np.diff(so.eqsys.W_cold[:,0]) / np.diff(t)
    sdWidt = np.diff(so.eqsys.W_i['D'][:,0]) / np.diff(t)
    sdIpdt = np.diff(so.eqsys.I_p[:,0]) / np.diff(t)
    sdIMK2dt = np.diff(so.eqsys.I_wall[:,0]) / np.diff(t)
    sdnD0dt = np.diff(so.eqsys.n_i['D'][0][:,0]) / np.diff(t)
    sdnD1dt = np.diff(so.eqsys.n_i['D'][1][:,0]) / np.diff(t)

    for i in range(1, t.size):
        x = fromSTREAM(so, uqh, time=i)

        dWedt.append(equations[uqh.map['We']](t[i], x)[0])
        dWidt.append(equations[uqh.map['Wi']](t[i], x)[0])
        dIpdt.append(equations[uqh.map['Ip']](t[i], x)[0])
        dIMK2dt.append(equations[uqh.map['IMK2']](t[i], x))
        dnD0dt.append(equations[uqh.map['niD_0']](t[i], x)[0])
        dnD1dt.append(equations[uqh.map['niD_1']](t[i], x)[0])

    # Plot results
    plotTimeDerivatives(t[1:], sdWedt, dWedt, title='Electron heat')
    plotTimeDerivatives(t[1:], sdWidt, dWidt, title='Ion heat')
    plotTimeDerivatives(t[1:], sdIpdt, dIpdt, title='Plasma current')
    plotTimeDerivatives(t[1:], sdIMK2dt, dIMK2dt, title='MK2 current')
    plotTimeDerivatives(t[1:], sdnD0dt, dnD0dt, title='D-0 density')
    plotTimeDerivatives(t[1:], sdnD1dt, dnD1dt, title='D-1 density')

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
        varNames = term['pydyon'].__init__.__code__.co_varnames
    else:
        varNames = term['pydyon'].__code__.co_varnames
    s = _getSubSettings(settings, varNames)

    if 'quantities' in varNames:
        s['quantities'] = unknowns
    if 'ions' in varNames:
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
            y[i] = term['eval'](obj, t[i], x, unknowns)
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
    for ion in so.eqsys.n_i.ions:
        dct[f'ni{ion.name}']=ion.data[time,:,0]

    return uqh.setvector(dct,so.eqsys.grid.t[time])


def loadSTREAMSettings(ss):
    """
    Convert a STREAMSettings object into a dict with settings
    for PYDYON.
    """
    # Load settings from STREAMSettings object
    settings = {}
    for s in SETTINGS.keys():
        settings[s] = SETTINGS[s](ss)

    return settings


def plotTimeDerivatives(t, dt1, dt2, title=None):
    """
    Plot two time derivatives (helper function for 'compareToSTREAMdt()')
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(t, dt1, 'k-', label='STREAM')
    ax.plot(t, dt2, 'r--', label='PYDYON')

    ax.set_xlabel('Time (s)')

    if title is not None:
        ax.set_title(title)

    ax.legend()

    return ax


