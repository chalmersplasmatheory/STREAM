# Representation of the mean-free path

from DREAM.Output.UnknownQuantity import UnknownQuantity


class MeanFreePath(UnknownQuantity):
    

    def __init__(self, name, data, attr, grid, output):
        """
        Constructor.
        """
        super().__init__(name, data, attr, grid, output)
        

    def __repr__(self):
        """
        Convert this object to an "official" string.
        """
        s = self.__str__() + "\n"
        if hasattr(self, "description") and hasattr(self, "description_eqn"):
            s += ":: {}\n::Evolved using: {}\n".format(self.description, self.description_eqn)
        s += self.dumps()

        return s


    def __str__(self):
        """
        Convert this object to a string.
        """
        return '({}) Mean-free path of size NI x NT x NR = {} x {} x {}'.format(self.name, self.data.shape[0], self.data.shape[1], self.data.shape[2])


    def __getitem__(self, index):
        """
        Direct access to data.
        """
        return self.data[index]

    
    def get(self, ion=None, r=None, t=None):
        """
        Returns data for the specified ion, or in the specified time
        interval or radial point. If none of the indices are given, returns
        the full evolution of the quantity.
        """
        sion = ion if ion is not None else slice(None)
        sr = r if r is not None else slice(None)
        st = t if t is not None else slice(None)

        return self.data[sion,sr,st]


    def dumps(self, ion=None, r=None, t=None):
        """
        Print the data in this quantity.
        """
        return self.get(ion=ion, r=r, t=t).__str__()


