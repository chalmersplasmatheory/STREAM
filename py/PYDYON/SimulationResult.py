# Encapsulation of PYDYON simulation result

import matplotlib.pyplot as plt
import numpy as np


class SimulationResult:
    

    def __init__(self, t, x):
        """
        Constructor.
        """
        self.t = t
        self.x = x


    def plot(self):
        """
        Plot a summary of the results.
        """
        fig, axs = plt.subplots(2, 3, figsize=(10, 6))

        self.plotQuantity('ne', axs[0,0], color='k', title='Electron density')
        self.plotQuantity('Te', axs[0,1], color='k', title='Electron temperature')
        self.plotQuantity('Ti', axs[0,2], color='k', title='Ion temperature')

        nD = self.x['niD']
        self.plotQuantity(nD[0], axs[1,0], color='k', label='D-0', title='Deuterium')
        self.plotQuantity(nD[1], axs[1,0], color='r', label='D-1')
        axs[1,0].set_ylabel('nD')
        axs[1,0].legend()

        self.plotQuantity('Ip', axs[1,1], color='k', title='Plasma current')
        #self.plotQuantity('IMK2', axs[1,2], color='k', title='MK2 current')

        gamma_i = self.x['ne'] / (nD[0] + self.x['ne'])
        self.plotQuantity(gamma_i*100, axs[1,2], color='k', title='Ionization fraction')
        axs[1,2].set_ylabel(r'\%')
        axs[1,2].set_ylim([0, 105])

        plt.tight_layout()
        plt.show()


    def plotQuantity(self, data, ax=None, title=None, *args, **kwargs):
        """
        Plot a named unknown quantity.
        """
        if ax is None:
            ax = plt.axes()

        if type(data) == str:
            name = data
            data = self.x[data]
        else:
            name = None

        ax.plot(self.t, data, *args, **kwargs)
        ax.set_xlabel('Time (s)')
        if name is not None:
            ax.set_ylabel(f'{name}')

        if title is not None:
            ax.set_title(title)

        return ax


