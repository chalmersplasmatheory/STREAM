import numpy as np
#from .. DREAMException import DREAMException
from DREAM.Settings import TransportSettings as DREAMTransport

TRANSPORT_NONE = DREAMTransport.TRANSPORT_NONE
TRANSPORT_PRESCRIBED = DREAMTransport.TRANSPORT_PRESCRIBED
TRANSPORT_RECHESTER_ROSENBLUTH = DREAMTransport.TRANSPORT_RECHESTER_ROSENBLUTH
TRANSPORT_SVENSSON = DREAMTransport.TRANSPORT_SVENSSON
TRANSPORT_DYON = 5

INTERP3D_NEAREST = DREAMTransport.INTERP3D_NEAREST
INTERP3D_LINEAR  = DREAMTransport.INTERP3D_LINEAR

INTERP1D_NEAREST = DREAMTransport.INTERP1D_NEAREST
INTERP1D_LINEAR  = DREAMTransport.INTERP1D_LINEAR

SVENSSON_INTERP1D_PARAM_TIME = DREAMTransport.SVENSSON_INTERP1D_PARAM_TIME
SVENSSON_INTERP1D_PARAM_IP   = DREAMTransport.SVENSSON_INTERP1D_PARAM_IP

BC_CONSERVATIVE = DREAMTransport.BC_CONSERVATIVE     # Assume no flux through r=rmax
BC_F_0 = DREAMTransport.BC_F_0                       # Assume f=0 outside the plasma
BC_DF_CONST = DREAMTransport.BC_DF_CONST             # Assume that df/dr is constant on the plasma boundary

class TransportSettings(DREAMTransport.TransportSettings):
    def __init__(self, kinetic=False):
        super().__init__(kinetic=kinetic)
        
    def setDYON(self, t=None):
        self.type = TRANSPORT_DYON
        
    def verifySettings(self):
        if self.type == TRANSPORT_DYON:
            pass
        else:
            super().verifySettings()
