#
# Representation of the output from a STREAM simulation.
#########################################################


import numpy as np
from DREAM import DREAMOutput

class STREAMOutput(DREAMOutput):
    

    def __init__(self, filename=None, path="", lazy=True):
        """
        Construct a new ``STREAMOutput`` object. If ``filename`` is given, the
        object is read from the (HDF5) file with that name. If ``path`` is also
        given, this is used to locate the group in the file which contains the
        output.

        :param str filename: Name of file to load output from.
        :param str path:     Path to group in HDF5 file containing the output.
        :param bool lazy:    If ``True``, only loads data from the file when explicitly requested.
        """
        super().__init__(filename=filename, path=path, lazy=lazy)

