# Specialized ion charge state class for STREAM

import matplotlib.pyplot as plt
from DREAM.Output.IonState import IonState as DREAMIonState


class IonState(DREAMIonState):
    

    def __init__(self, name, Z, Z0, data, grid, output, attr=list()):
        """
        Constructor.
        """
        super().__init__(name=name, Z=Z, Z0=Z0, data=data, grid=grid, output=output, attr=attr)


    def plotBalance(self, ax=None, show=None):
        """
        Generates a plot showing the impact of the various terms in the
        ion rate equation (if saved).
        """
        if 'stream' not in self.output.other:
            raise Exception("The 'stream' category of other quantities was not saved to this output.")

        genax = ax is None

        if genax:
            ax = plt.axes()

            if show is None:
                show = True

        qty = {
            'ionrateequation_posIonization': 'Positive ionization',
            'ionrateequation_negIonization': 'Negative ionization',
            'ionrateequation_posRecombination': 'Positive recombination',
            'ionrateequation_negRecombination': 'Negative recombination',
            'ionrateequation_posChargeExchange': 'Positive charge-exchange',
            'ionrateequation_negChargeExchange': 'Negative charge-exchange'
        }

        nm = self.name.partition('-')[0]
        offset = self.output.eqsys.n_i.getIonOffset(nm, self.Z0)

        stream = self.output.other.stream
        for term, desc in qty.items():
            trm = stream[term][:,offset,0]
            trm = trm.reshape((trm.size,))
            ax.plot(self.grid.t[1:], trm, label=desc)

        ax.legend()

        if show:
            plt.show(block=False)

        return ax
            

