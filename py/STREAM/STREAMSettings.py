#
# An object representing the settings passed when running STREAM.
###################################################################

import copy
import numpy as np
from DREAM import DREAMIO, DREAMSettings

from .Settings.EquationSystem import EquationSystem
from .Settings.RadialGrid import RadialGrid


class STREAMSettings(DREAMSettings):
    

    def __init__(self, filename=None, path="", chain=True, keepignore=False):
        """
        Construct a new STREAMSettings object. If ``filename`` is given,
        the object is read from the (HDF5) file with that name.
        If ``path`` is also given, this is used to locate the group
        in the file which contains the settings.

        :param str filename:    Name of the file to load settings from.
        :param str path:        Path to group in HDF5 file containing the settings.
        :param bool chain:      If ``True``, sets the newly created ``STREAMSettings`` object to take the output of the simulation contained in the STREAM output file ``filename`` as input (i.e. calls :py:method:`fromOutput`).
        :param bool keepignore: If ``True``, keeps the list of unknown quantities to ignore when initializing from a previous simulation (and ``chain=True``).
        """
        super().__init__(filename=None)

        # TODO replace 'EquationSystem'
        self.addSetting('radialgrid', RadialGrid())
        
        # Should be defined last as it may need access to the
        # obejcts created above...
        self.addSetting('eqsys', EquationSystem(settings=self))

        if filename is not None:
            if type(filename) == str:
                self.load(filename, path=path, lazy=False)
            elif type(filename) == STREAMSettings:
                td = filename.todict()
                self.fromdict(td)

                if chain:
                    self.fromOutput(filename.output.filename)
                    self.output.setFilename('output.h5')

                    if not keepignore:
                        self.clearIgnore()
