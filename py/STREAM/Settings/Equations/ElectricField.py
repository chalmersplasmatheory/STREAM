# Electric field settings for STREAM

from DREAM.Settings.Equations.ElectricField import ElectricField as DREAMEfield
from DREAM.Settings.Equations.PrescribedScalarParameter import PrescribedScalarParameter


TYPE_PRESCRIBED = 1
TYPE_SELFCONSISTENT = 2
TYPE_CIRCUIT = 3


class ElectricField(DREAMEfield):
    

    def __init__(self, ttype=TYPE_CIRCUIT, *args, **kwargs):
        """
        Constructor.
        """
        super().__init__(*args, **kwargs)
        
        self.setType(ttype)

        self.circuit_Lp = 0
        self.circuit_Lwall = 0
        self.circuit_M = 0
        self.circuit_Rwall = 0

        self.circuit_Vloop = None
        self.circuit_Vloop_t = None


    def setType(self, ttype):
        if ttype == TYPE_CIRCUIT:
            self.type = ttype
        else:
            super().setType(ttype)


    def setInductances(self, Lp, Lwall, M, Rwall):
        """
        Set the self and mutual inductances of the plasma and wall, as well
        as the resistance of the wall. Used only for the electric field type
        ``TYPE_CIRCUIT``.
        """
        if Lp <= 0:
            raise Exception("Invalid value assigned to 'Lp'.")
        elif Lwall <= 0:
            raise Exception("Invalid value assigned to 'Lwall'.")
        elif M <= 0:
            raise Exception("Invalid value assigned to 'M'.")
        elif Rwall <= 0:
            raise Exception("Invalid value assigned to 'Rwall'.")

        self.circuit_Lp = Lp
        self.circuit_Lwall = Lwall
        self.circuit_M = M
        self.circuit_Rwall = Rwall


    def setCircuitVloop(self, Vloop, times=0):
        """
        When ``TYPE_CIRCUIT``, evolves the external loop voltage according
        to this value.
        """
        _data, _time = self._setScalarData(data=Vloop, times=times)
        self.circuit_Vloop = _data
        self.circuit_Vloop_t = _time


    def fromdict(self, data):
        """
        Sets this parameter from settings provided in the given dictionary.
        """
        self.type = data['type']

        if self.type == TYPE_CIRCUIT:
            self.efield = data['init']['x']
            self.radius = data['init']['r']

            self.circuit_Lp = data['circuit']['Lp']
            self.circuit_Lwall = data['circuit']['Lwall']
            self.circuit_M = data['circuit']['M']
            self.circuit_Rwall = data['circuit']['Rwall']

            self.circuit_Vloop = data['circuit']['Vloop']['x']
            self.circuit_Vloop_t = data['circuit']['Vloop']['t']

            self.verifySettings()
        else:
            super().fromdict(data)


    def todict(self):
        """
        Returns a Python dictionary containing all settings of this
        ElectricField object.
        """
        if self.type == TYPE_CIRCUIT:
            data = {
                'type': self.type,
                'init': {
                    'x': self.efield,
                    'r': self.radius
                },
                'circuit': {
                    'Lp': self.circuit_Lp,
                    'Lwall': self.circuit_Lwall,
                    'M': self.circuit_M,
                    'Rwall': self.circuit_Rwall,
                    'Vloop': {
                        'x': self.circuit_Vloop,
                        't': self.circuit_Vloop_t
                    }
                }
            }

            return data
        else:
            return super().todict()

    
    def verifySettings(self):
        """
        Verify that the settings of this unknown are correctly set.
        """
        if self.type == TYPE_CIRCUIT:
            if self.circuit_Lp <= 0:
                raise Exception("E_field: Parameter 'Lp' must be given a value > 0.")
            elif self.circuit_Lwall <= 0:
                raise Exception("E_field: Parameter 'Lwall' must be given a value > 0.")
            elif self.circuit_M <= 0:
                raise Exception("E_field: Parameter 'M' must be given a value > 0.")
            elif self.circuit_Rwall <= 0:
                raise Exception("E_field: Parameter 'Rwall' must be given a value > 0.")

            PrescribedScalarParameter._verifySettingsPrescribedScalarData(self, name='Vloop', data=self.circuit_Vloop, times=self.circuit_Vloop_t)
            self._verifySettingsPrescribedInitialData()
        else:
            super().verifySettings()


