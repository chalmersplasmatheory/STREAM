
from DREAM.Settings.Equations.PrescribedParameter import PrescribedParameter

class RadialGrid(PrescribedParameter):
    

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
        self.vessel_volume = None
        self.c1 = None
        self.c2 = None
        self.c3 = None


    def setB0(self, B0, t=0):
        """
        Prescribe the (time evolution of the) magnetic field strength
        on the magnetic axis.

        :param B0: Magnetic field strength on the magnetic axis.
        :param t:  Time vector (if ``B0`` varies with time).
        """
        self.B0, _, self.tB0 = self._setPrescribedData(data=B0, times=t)


    def setMinorRadius(self, a, t=0):
        """
        Prescribe the (time evolution of the) plasma minor radius.

        :param a: Plasma minor radius.
        :param t: Time vector (if ``a`` varies with time).
        """
        self.a, _, self.ta = self._setPrescribedData(data=a, times=t)


    def setElongation(self, kappa, t=0):
        """
        Prescribe the (time evolution of the) plasma elongation parameter.

        :param kappa: Plasma elongation.
        :param t:     Time vector (if ``kappa`` varies with time).
        """
        self.kappa, _, self.tkappa = self._setPrescribedData(data=kappa, times=t)
        
    def setVesselVolume(self, v):
        """
        Prescribe the vacuum vessel volume.

        :param v: Vacuum vessel volume.
        """
        self.vessel_volume = v

    def setRecyclingCoefficient1(self, c1):
        """
        Prescribe the first recycling coefficient for deuterium.

        :param c1: recycling coefficient.
        """
        self.c1 = c1
		
    def setRecyclingCoefficient2(self, c2):
        """
        Prescribe the second recycling coefficient for deuterium.

        :param c2: recycling coefficient.
        """
        self.c2 = c2
		
    def setRecyclingCoefficient3(self, c3):
        """
        Prescribe the third recycling coefficient for deuterium.

        :param c3: recycling coefficient.
        """
        self.c3 = c3
	
    def fromdict(self, data):
        """
        Load settings from the given dictionary.
        """
        self.a, self.ta = data['a']['x'], data['a']['t']
        self.B0, self.tB0 = data['B0']['x'], data['B0']['t']
        self.kappa, self.tkappa = data['kappa']['x'], data['kappa']['t']
        self.vessel_volume = data['vessel_volume']['x']
        self.c1 = data['c1']
        self.c2 = data['c2']
        self.c3 = data['c3']


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
            'vessel_volume': {
                'x': self.vessel_volume
            },
            'c1': self.c1,
            'c2': self.c2,
            'c3': self.c3
        }

        return data


    def verifySettings(self):
        """
        Verify that the RadialGrid settings are consistent.
        """
        r0 = np.asarray([0])

        self._verifySettingsPrescribedData('a', self.a, r0, self.ta)
        self._verifySettingsPrescribedData('B0', self.B0, r0, self.tB0)
        self._verifySettingsPrescribedData('kappa', self.kappa, r0, self.tkappa)
        if type(self.vessel_volume) != float:
            raise TypeError('The prescribed vessel volume must be of type float')
        if type(self.c1) != float:
            raise TypeError('The prescribed recycle coefficient 1 must be of type float')
        if type(self.c2) != float:
            raise TypeError('The prescribed recycle coefficient 2 must be of type float')
        if type(self.c3) != float:
            raise TypeError('The prescribed recycle coefficient 3 must be of type float')

