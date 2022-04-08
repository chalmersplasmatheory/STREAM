import numpy as np

from DREAM.Settings.Equations.PrescribedScalarParameter import PrescribedScalarParameter

class RadialGrid(PrescribedScalarParameter):
    

    def __init__(self):
        """
        Constructor.
        """
        # Only here for compatibility with DREAM parts of
        # the settings interface
        self.nr = 1

        self.a, self.ta = None, None
        self.B0, self.tB0 = None, None
        self.kappa, self.tkappa = None, None
        self.delta, self.tdelta = None, None
        self.vessel_volume = None
        self.b = 0.0
        self.R0 = 2.0
        self.Iref = 100.0e3
        self.Bv = 1.0e-3
        self.connectionLengthFactor = 3.0
        self.P_inj = 0.0
        self.f_o = 0.5
        self.f_x = 0.5
        self.theta = np.pi/4
        self.phi = 0.0
        self.N = 1

        self.setElongation(1)
        self.setTriangularity(0)



    def setB0(self, B0, t=0):
        """
        Prescribe the (time evolution of the) magnetic field strength
        on the magnetic axis.

        :param B0: Magnetic field strength on the magnetic axis.
        :param t:  Time vector (if ``B0`` varies with time).
        """
        self.B0, self.tB0 = self._setScalarData(data=B0, times=t)


    def setBv(self, Bv):
        """
        Sets the value of the stray magnetic field.
        """
        self.Bv = float(Bv)

    def setConnectionLengthFactor(self, connectionLengthFactor):
        """
        Sets the value of the stray magnetic field.
        """
        self.connectionLengthFactor = float(connectionLengthFactor)


    def setIref(self, Iref):
        """
        Set the value of the reference plasma current at which
        closed flux surfaces are assumed to form.
        """
        self.Iref = float(Iref)


    def setMinorRadius(self, a, t=0):
        """
        Prescribe the (time evolution of the) plasma minor radius.

        :param a: Plasma minor radius.
        :param t: Time vector (if ``a`` varies with time).
        """
        self.a, self.ta = self._setScalarData(data=a, times=t)


    def setMajorRadius(self, R0):
        """
        (Analytic toroidal)
        Set the tokamak major radius.
        """
        if R0 <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to major radius 'R0': {}".format(R0))

        self.R0 = float(R0)


    def setWallRadius(self, wall_radius):
        """
        (Cylindrical, Analytic toroidal)
        Set the minor radius of the wall
        """
        self.b = float(wall_radius)


    def setElongation(self, kappa, t=0):
        """
        Prescribe the (time evolution of the) plasma elongation parameter.

        :param kappa: Plasma elongation.
        :param t:     Time vector (if ``kappa`` varies with time).
        """
        self.kappa, self.tkappa = self._setScalarData(data=kappa, times=t)

    def setTriangularity(self, delta, t=0):
        """
        Prescribe the (time evolution of the) plasma triangularity parameter.

        :param delta: Plasma elongation.
        :param t:     Time vector (if ``delt`` varies with time).
        """
        self.delta, self.tdelta = self._setScalarData(data=delta, times=t)

    def setVesselVolume(self, v):
        """
        Prescribe the vacuum vessel volume.

        :param v: Vacuum vessel volume.
        """
        self.vessel_volume = float(v)

    def setInjectedECHPower(self, P_inj):
        """
        Prescribe the injected ECH power.

        :param P_inj: injected ECH power.
        """
        self.P_inj = float(P_inj)

    def setFractionOMode(self, f_o):
        """
        Prescribe the fraction of ECH that is O mode.

        :param f_o: fraction of O mode.
        """
        self.f_o = float(f_o)

    def setFractionXMode(self, f_x):
        """
        Prescribe the fraction of ECH that is X mode.

        :param f_x: fraction of X mode.
        """
        self.f_x = float(f_x)

    def setECHPoloidalAngle(self, theta):
        """
        Prescribe poloidal angle between EC-beam path and vertical z-axis.

        :param theta: poloidal angle.
        """
        self.theta = float(theta)

    def setECHPoloidalAngle(self, phi):
        """
        Prescribe toroidal angle between magnetic field and EC wave injection.

        :param phi: toroidal angle.
        """
        self.phi = float(phi)

    def setFundamentalHarmonic(self, N):
        """
        Prescribe fundamental harmonic of ECH beam.

        :param theta: fundamental harmonic.
        """
        self.N = int(N)

    def setECHParameters(self, P_inj, f_o, f_x, theta, phi, N):
        """
        Prescribe the injected ECH power.

        :param P_inj: injected ECH power.
        :param f_o: fraction of O mode.
        :param f_x: fraction of X mode.
        :param theta: poloidal angle.
        :param phi: toroidal angle.
        :param theta: fundamental harmonic.
        """
        self.P_inj = float(P_inj)
        self.f_o   = float(f_o)
        self.f_x   = float(f_x)
        self.theta = float(theta)
        self.phi   = float(phi)
        self.N     = int(N)

    def fromdict(self, data):
        """
        Load settings from the given dictionary.
        """
        self.a, self.ta = data['a']['x'], data['a']['t']
        self.B0, self.tB0 = data['B0']['x'], data['B0']['t']
        self.kappa, self.tkappa = data['kappa']['x'], data['kappa']['t']
        self.delta, self.tdelta = data['delta']['x'], data['delta']['t']
        self.vessel_volume = data['wall']['vessel_volume']
        self.R0 = data['R0']

        if 'wall_radius' in data:
            self.b = data['wall_radius']
            if type(self.b) == np.ndarray:
                self.b = float(self.b[0])
            else:
                self.b = float(self.b)

        if 'Bv' in data:
            self.Bv = data['Bv']

        if 'connectionLengthFactor' in data:
            self.connectionLengthFactor = data['connectionLengthFactor']

        if 'Iref' in data:
            self.Iref = data['Iref']

        if 'P_inj' in data:
            self.P_inj = data['P_inj']

        if 'f_o' in data:
            self.f_o = data['f_o']

        if 'f_x' in data:
            self.f_x = data['f_x']

        if 'theta' in data:
            self.theta = data['theta']

        if 'phi' in data:
            self.phi = data['phi']

        if 'N' in data:
            self.N = data['N']



    def todict(self, verify=True):
        """
        Returns the settings in this object as a Python dictionary.
        """
        if verify:
            self.verifySettings()

        data = {
            'a': {
                't': self.ta,
                'x': self.a
            },
            'B0': {
                't': self.tB0,
                'x': self.B0
            },
            'kappa': {
                't': self.tkappa,
                'x': self.kappa
            },
            'delta': {
                't': self.tdelta,
                'x': self.delta
            },
            'Bv': self.Bv,
            'connectionLengthFactor': self.connectionLengthFactor,
            'Iref': self.Iref,
            'R0': self.R0,
            'wall_radius': self.b,
            'wall' : {
                'vessel_volume': self.vessel_volume
            },
            'P_inj': self.P_inj,
            'f_o': self.f_o,
            'f_x': self.f_x,
            'theta': self.theta,
            'phi': self.phi,
            'N': self.N
        }

        return data


    def verifySettings(self):
        """
        Verify that the RadialGrid settings are consistent.
        """
        self._verifySettingsPrescribedScalarData('a', self.a, self.ta)
        self._verifySettingsPrescribedScalarData('B0', self.B0, self.tB0)
        self._verifySettingsPrescribedScalarData('kappa', self.kappa, self.tkappa)
        self._verifySettingsPrescribedScalarData('delta', self.delta, self.tdelta)

        if type(self.vessel_volume) != float:
            raise TypeError('The prescribed vessel volume must be of type float')
        if type(self.Bv) != float:
            raise TypeError('The prescribed stray magnetic field must be of type float')
        if type(self.connectionLengthFactor) != float:
            raise TypeError('The prescribed connection length factor must be of type float')
        if type(self.Iref) != float:
            raise TypeError('The prescribed reference current must be of type float')
        if type(self.Iref) != float:
            raise TypeError('The prescribed reference plasma current must be of type float')
        if type(self.P_inj) != float:
            raise TypeError('The prescribed injected ECH power must be of type float')
        if type(self.f_o) != float:
            raise TypeError('The fraction of O mode must be of type float')
        if type(self.f_x) != float:
            raise TypeError('The fraction of X mode must be of type float')
        if type(self.theta) != float:
            raise TypeError('The poloidal angle of ECH-beam must be of type float')
        if type(self.phi) != float:
            raise TypeError('The toroidal angle of ECH-wave injection must be of type float')
        if type(self.N) != int:
            raise TypeError('The fundmental harmonic of ECH beam must be of type int')


        if self.R0 is None or self.R0 <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to tokamak major radius 'R0': {}".format(self.R0))
        if not np.isscalar(self.b):
            raise DREAMException("RadialGrid: The specified wall radius is not a scalar: {}.".format(self.b))


        if self.f_o + self.f_x != 1 and self.f_o + self.f_x != 0:
            self.f_o = self.f_o / (self.f_o + self.f_x)
            self.f_x = self.f_x / (self.f_o + self.f_x)
            print('WARNING: fractions of O and X mode modified to sum to 1.')
        if (self.theta - np.pi/2) % np.pi == 0 or self.theta > 2*np.pi or self.theta <= 0:
            raise DREAMException(
                "RadialGrid: Invalid value assigned to poloidal angle of ECH-beam 'theta': {}".format(self.theta))


