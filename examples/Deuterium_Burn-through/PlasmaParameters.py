
import numpy as np
import scipy.constants
import scipy.interpolate
import scipy.special


def evaluateSpitzerConductivity(n, T, Z):
    """
    Evaluate the Spitzer conductivity at the given temperature and charge.
    Uses the formula on p. 72 in Helander & Sigmar.
    """
    sigma = getSpitzerConductivity()
    return sigma(n, T, Z)

def getCoulombLogarithm(T, n):
    """
    Calculates the Coulomb logarithm according to the formula given in
    Wesson's book "Tokamaks".

    :param float T: Plasma temperature (eV).
    :param float n: Plasma density (m^-3).
    """
    return 14.9 - 0.5*np.log(n / 1e20) + np.log(T / 1e3)


def getEc(T, n):
    """
    Calculates the Connor-Hastie critical electric field, below which no
    runaway electrons can be generated.

    :param float T: Plasma temperature (eV).
    :param float n: Plasma density (m^-3).
    """
    logLambda = getCoulombLogarithm(T, n)

    c = scipy.constants.c
    e = scipy.constants.e
    me = scipy.constants.m_e
    eps0 = scipy.constants.epsilon_0

    return (n*logLambda*e**3) / (4*np.pi*eps0**2 * me * c**2)


def getConnorHastieCriticalField(T, n): return getEc(T, n)


def getED(T, n):
    """
    Calculates the Dreicer electric field at the given plasma temperature and
    density, giving the electric field at which all electrons are accelerated
    to the runaway region.

    :param float T: Plasma temperature (eV).
    :param float n: Plasma density (m^-3).
    """
    c = scipy.constants.c
    e = scipy.constants.e
    me = scipy.constants.m_e

    Ec = getEc(T, n)

    return Ec * me * c**2 / (e*T)


def getDreicerElectricField(T, n): return getED(T, N)


def getSpitzerConductivity():
    """
    Returns a function corresponding to the Spitzer conductivity.
    Uses the formula on p. 72 in Helander & Sigmar.
    """
    eps0 = scipy.constants.epsilon_0
    e    = scipy.constants.e
    me   = scipy.constants.m_e

    # The prefactor has a Z dependence which is only known numerically
    # (table on p. 74 in Helander & Sigmar and Table III in Spitzer & Härm).
    # Therefore we interpolate to find the value for any Z. The data
    # contains a point at Z=inf, so we interpolate in 1/Z instead.
    invZs  = np.array([0, 1/16, 0.25, 0.5, 1])
    vals   = 32/(3*np.pi)*np.array([1, 0.9225, 0.7849, 0.6833, 0.5816])
    preFac = scipy.interpolate.PchipInterpolator(invZs, vals)
    
    pf = 3/(4*np.sqrt(2*np.pi)*e**2*np.sqrt(me)) * (4*np.pi*eps0)**2

    return lambda n, T, Z : preFac(1/Z) * pf * (T*e)**(3/2) / (Z*getCoulombLogarithm(T, n))


def getTauEETh(T, n):
    """
    Calculates the thermal electron-electron collision frequency.

    :param float T: Plasma temperature (eV).
    :param float n: Plasma density (m^-3).
    """
    mc2  = scipy.constants.physical_constants['electron mass energy equivalent in MeV'][0] * 1e6
    betaTh2 = getBetaThermal(T)**2

    return getTauEERel(T, n) * betaTh2*np.sqrt(betaTh2)


def getThermalElectronCollisionFrequency(T, n): return getTauEETh(T, n)


def getBetaThermal(T): return getNormalizedThermalSpeed(T)


def getNormalizedThermalSpeed(T):
    """
    Calculates the thermal electron speed, v_{th}, normalized
    to the speed of light in vacuum.
    """
    return getThermalSpeed(T) / scipy.constants.c


def getThermalSpeed(T):
    """
    Calculates the thermal electron speed, v_{th}, for the given
    plasma temperature.

    :param float T: Electron temperature.
    """
    return np.sqrt(2*scipy.constants.e*T / scipy.constants.m_e)


def getTauEERel(T, n):
    """
    Calculates the relativistic electron-electron collision frequency.

    :param float T: Plasma temperature (eV).
    :param float n: Plasma density (m^-3).
    """
    c = scipy.constants.c
    e = scipy.constants.e
    me = scipy.constants.m_e
    r0 = scipy.constants.physical_constants['classical electron radius'][0]

    C = 4*np.pi * r0**2 * c
    logLambda = getCoulombLogarithm(T, n)

    return 1/(logLambda * n * C)


def getRelativisticElectronCollisionFrequency(T, n): return getTauEERel(T, n)


