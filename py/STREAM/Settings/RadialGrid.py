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
        self.c1 = None
        self.c2 = None
        self.c3 = None
        self.b = 0.0
        self.R0 = 2.0
        self.Iref = 100.0e3
        self.Bv = 1.0e-3
        self.uncertaintyFactor = 1.0

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


    def setRunawayConfinementUncertainty(self, uncertaintyFactor):
        """
        Sets the value of the runaway confinement factor
        """
        self.uncertaintyFactor = float(uncertaintyFactor)



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

    def setRecyclingCoefficient1(self, c1):
        """
        Prescribe the first recycling coefficient for deuterium.

        :param c1: recycling coefficient.
        """
        self.c1 = float(c1)
		
    def setRecyclingCoefficient2(self, c2):
        """
        Prescribe the second recycling coefficient for deuterium.

        :param c2: recycling coefficient.
        """
        self.c2 = float(c2)
		
    def setRecyclingCoefficient3(self, c3):
        """
        Prescribe the third recycling coefficient for deuterium.

        :param c3: recycling coefficient.
        """
        self.c3 = float(c3)
	
    def fromdict(self, data):
        """
        Load settings from the given dictionary.
        """
        self.a, self.ta = data['a']['x'], data['a']['t']
        self.B0, self.tB0 = data['B0']['x'], data['B0']['t']
        self.kappa, self.tkappa = data['kappa']['x'], data['kappa']['t']
        self.delta, self.tdelta = data['delta']['x'], data['delta']['t']
        self.vessel_volume = data['wall']['vessel_volume']
        self.c1 = data['wall']['c1']
        self.c2 = data['wall']['c2']
        self.c3 = data['wall']['c3']
        self.R0 = data['R0']

        if 'wall_radius' in data:
            self.b = data['wall_radius']
            if type(self.b) == np.ndarray:
                self.b = float(self.b[0])
            else:
                self.b = float(self.b)
        
        if 'Bv' in data:
            self.Bv = data['Bv']

        if 'uncertaintyFactor' in data:
            self.uncertaintyFactor = data['uncertaintyFactor']


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
            'uncertaintyFactor': self.uncertaintyFactor,
            'R0': self.R0,
            'wall_radius': self.b,
            'wall' : {
                'vessel_volume': self.vessel_volume,
                'c1': self.c1,
                'c2': self.c2,
                'c3': self.c3
            }
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
        if type(self.c1) != float:
            raise TypeError('The prescribed recycle coefficient 1 must be of type float')
        if type(self.c2) != float:
            raise TypeError('The prescribed recycle coefficient 2 must be of type float')
        if type(self.c3) != float:
            raise TypeError('The prescribed recycle coefficient 3 must be of type float')
        if type(self.Bv) != float:
            raise TypeError('The prescribed stray magnetic field must be of type float')
        if type(self.uncertaintyFactor) != float:
            raise TypeError('The prescribed runaway confinement uncertainty factor must be of type float')
        if type(self.Iref) != float:
            raise TypeError('The prescribed reference plasma current must be of type float')
            
        if self.R0 is None or self.R0 <= 0:
            raise DREAMException("RadialGrid: Invalid value assigned to tokamak major radius 'R0': {}".format(self.R0))
        if not np.isscalar(self.b):
            raise DREAMException("RadialGrid: The specified wall radius is not a scalar: {}.".format(self.b))


