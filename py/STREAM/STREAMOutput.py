#
# Representation of the output from a STREAM simulation.
#########################################################


import numpy as np
from DREAM import DREAMIO, DREAMOutput

from DREAM.Output.OtherIonSpeciesFluidQuantity import OtherIonSpeciesFluidQuantity
from DREAM.Output.OtherIonSpeciesScalarQuantity import OtherIonSpeciesScalarQuantity
from . Output.IonHandler import IonHandler
from . Output.MeanFreePath import MeanFreePath
from DREAM.Output.ScalarQuantity import ScalarQuantity

from . STREAMSettings import STREAMSettings


class STREAMOutput(DREAMOutput):
    

    def __init__(self, filename=None, path="", lazy=True, loadsettings=True):
        """
        Construct a new ``STREAMOutput`` object. If ``filename`` is given, the
        object is read from the (HDF5) file with that name. If ``path`` is also
        given, this is used to locate the group in the file which contains the
        output.

        :param str filename: Name of file to load output from.
        :param str path:     Path to group in HDF5 file containing the output.
        :param bool lazy:    If ``True``, only loads data from the file when explicitly requested.
        """
        super().__init__(filename=filename, path=path, lazy=lazy, loadsettings=loadsettings)


    def load(self, filename, path="", lazy=True, loadsettings=True, *args, **kwargs):
        """
        Loads STREAM output from the specified file. If 'path' is
        given, this indicates which group path in the file to load
        the output from.

        :param str filename: Name of file to load output from.
        :param str path:     Path to subset of HDF5 file containing STREAM output.
        :param bool lazy:    If ``True``, allows the file to be read lazily (on-demand) by returning h5py DataSet objects instead of the actual data (wrapped in a DREAM.DataObject).
        """
        od = super().load(filename=filename, path=path, lazy=lazy, loadsettings=False, *args, **kwargs)

        if 'settings' in od and loadsettings:
            s = od['settings']
            if lazy:
                s = DREAMIO.unlazy(s)

            self.settings = STREAMSettings(s)

        self.eqsys.resetUnknown('lambda_i', MeanFreePath)
        self.eqsys.resetUnknown('n_i', IonHandler)
        self.eqsys.resetUnknown('Vloop', ScalarQuantity)

        if self.other is not None and 'stream' in self.other:
            self.other.stream.resetQuantity('V_n', OtherIonSpeciesScalarQuantity)
            self.other.stream.resetQuantity('V_n_tot', OtherIonSpeciesScalarQuantity)
            self.other.stream.resetQuantity('neutralinflux', OtherIonSpeciesScalarQuantity)
            #self.other.stream.resetQuantity('iontransport', OtherIonSpeciesScalarQuantity)
        
