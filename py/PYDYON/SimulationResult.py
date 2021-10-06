# Encapsulation of PYDYON simulation result

import matplotlib.pyplot as plt
import numpy as np
from . Equations import *


class SimulationResult:
    

    def __init__(self, t, x, simulation):
        """
        Constructor.
        """
        self.t = t
        self.x = x

        self.simulation = simulation

    def getT(self):
        return self.t

    def evaluateTerm(self, term):
        """
        Evaluates a term of the given type.
        """
        uqh = self.simulation.unknowns
        val = []
        for i in range(self.t.size):
            dct = {
                'We': self.x['We'][i],
                'Wi': self.x['Wi'][i],
                'Ip': self.x['Ip'][i],
                'IMK2': self.x['IMK2'][i],
                'niD': self.x['niD'][:,i],
                'niO': self.x['niO'][:,i],
                'niC': self.x['niC'][:,i]
            }
            x = uqh.setvector(dct,t=self.t[i])
            val.append(term(self.t[i], x))

        return np.array(val)

    def evaluateTermIonState(self, term, ionname, Z0):
        """
        Evaluates a term of the given type.
        """
        uqh = self.simulation.unknowns
        val = []
        for i in range(self.t.size):
            dct = {
                'We': self.x['We'][i],
                'Wi': self.x['Wi'][i],
                'Ip': self.x['Ip'][i],
                'IMK2': self.x['IMK2'][i],
                'niD': self.x['niD'][:,i],
                'niO': self.x['niO'][:,i],
                'niC': self.x['niC'][:,i]
            }
            x = uqh.setvector(dct,t=self.t[i])
            val.append(term(self.t[i], x, ionname, Z0))

        return np.array(val)

    def evaluateTermFull(self, term, eq,  full=False):
        """
        Evaluates a term of the given type.
        """
        uqh = self.simulation.unknowns
        val = []
        for i in range(self.t.size):
            dct = {
                'We': self.x['We'][i],
                'Wi': self.x['Wi'][i],
                'Ip': self.x['Ip'][i],
                'IMK2': self.x['IMK2'][i],
                'niD': self.x['niD'][:,i],
                'niO': self.x['niO'][:,i],
                'niC': self.x['niC'][:,i]
            }
            x = uqh.setvector(dct,t=self.t[i])
            val.append(term(self.t[i], x, full)[eq])

        return np.array(val)


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

    def plotJET(self):
        """
        Plot a summary of the results.
        """
        fig, axs = plt.subplots(3, 2, figsize=(10, 6))

        self.plotQuantity('Ip', axs[0,0], color='k', title='Plasma current', axis=[0,0.3,0,6e5])
        self.plotQuantity('Te', axs[1, 0], color='k', title='Electron temperature', axis=[0,0.3,0,400])
        self.plotQuantity('Ti', axs[2, 0], color='k', title='Ion temperature', axis=[0,0.3,0,400])
        self.plotQuantity('ne', axs[0,1], color='k', title='Electron density', axis=[0,0.3,0,6e18])

        #nD = self.x['niD']
        #self.plotQuantity(nD[0], axs[1,0], color='k', label='D-0', title='Deuterium')
        #self.plotQuantity(nD[1], axs[1,0], color='r', label='D-1')
        #axs[1,0].set_ylabel('nD')
        #axs[1,0].legend()


        #self.plotQuantity('IMK2', axs[1,2], color='k', title='MK2 current')

        #gamma_i = self.x['ne'] / (nD[0] + self.x['ne'])
        #self.plotQuantity(gamma_i*100, axs[1,2], color='k', title='Ionization fraction')
        #axs[1,2].set_ylabel(r'\%')
        #axs[1,2].set_ylim([0, 105])

        plt.tight_layout()
        plt.show()


    def plotKimThesis45(self):
        """
        Create a plot resembling figure 4.5 Kim's PhD thesis.
        """
        fig, axs = plt.subplots(3, 1, figsize=(5, 8))
        Ip = self.x['Ip']

        self.plotQuantity(self.x['Ip']/1e4, axs[0], title='Plasma current')
        axs[0].set_xlim([0, self.t[-1]])
        axs[0].set_ylim([0, 10])

        nD = self.x['niD']
        gamma_i = self.x['ne'] / (nD[0] + self.x['ne'])
        self.plotQuantity(gamma_i*100, axs[1], title='Ionization fraction')
        axs[1].set_xlim([0, self.t[-1]])
        axs[1].set_ylim([0, 105])

        sim = self.simulation
        s = sim.settings
        tausettings = {'Bphi': s['Bphi'], 'Bv': s['Bv'], 'l_MK2': s['l_MK2']}

        Poh    = OhmicPowerTerm(sim.unknowns, sim.ions)
        Prad   = RadiatedPowerTerm(sim.unknowns, sim.ions)
        Pequi  = EquilibrationPowerTerm(sim.unknowns, sim.ions)
        Pconve = ElectronConvectivePowerTerm(sim.unknowns, sim.ions, **tausettings)

        Vp = sim.unknowns.getV_p()
        print('Plasma volume Vp = {} m^3'.format(Vp))
        vPoh = self.evaluateTerm(Poh) * Vp
        vPrad = self.evaluateTerm(Prad) * Vp
        vPequi = self.evaluateTerm(Pequi) * Vp
        vPconve = self.evaluateTerm(Pconve) * Vp
        vPtot = vPrad+vPequi+vPconve
        vPnet = vPoh - vPtot

        axs[2].plot(self.t, vPtot/1e6, 'r--', label='Total electron power loss')
        axs[2].plot(self.t, vPrad/1e6, 'b.-', label='Radiation+ionization')
        axs[2].plot(self.t, vPequi/1e6, 'g--', label='Equilibration')
        axs[2].plot(self.t, vPconve/1e6, 'm:', label='Electron transport')
        axs[2].plot(self.t, vPnet/1e6, 'k-', label='Net electron heating power')
        axs[2].set_title('Power balance')
        axs[2].legend(frameon=False)

        axs[2].set_xlim([0, self.t[-1]])
        axs[2].set_ylim([0, 0.5])

        fig.tight_layout()
        plt.show()


    def plotKim2020(self):
        fig, axs = plt.subplots(1, 2, figsize=(8, 5))

        sim = self.simulation
        s = sim.settings

        tausettings = {'Bphi': s['Bphi'], 'Bv': s['Bv'], 'l_MK2': s['l_MK2']}
        Poh    = OhmicPowerTerm(sim.unknowns, sim.ions)
        Prad   = RadiatedPowerTerm(sim.unknowns, sim.ions)
        Pequi  = EquilibrationPowerTerm(sim.unknowns, sim.ions)
        Pconve = ElectronConvectivePowerTerm(sim.unknowns, sim.ions, **tausettings)
        Pcx    = ChargeExchangePowerTerm(sim.unknowns, sim.ions)
        Pconvi = IonConvectivePowerTerm(sim.unknowns, sim.ions, **tausettings)

        Vp = sim.unknowns.getV_p()
        vPoh = self.evaluateTerm(Poh) * Vp
        vPrad = self.evaluateTerm(Prad) * Vp
        vPequi = self.evaluateTerm(Pequi) * Vp
        vPconve = self.evaluateTerm(Pconve) * Vp
        vPtot = vPrad+vPequi+vPconve
        vPnet = vPoh - vPtot

        vPcx = self.evaluateTerm(Pcx) * Vp
        vPconvi = self.evaluateTerm(Pconvi) * Vp

        #axs[0].plot(self.t, vPtot/1e6, 'r--', label='Total electron power loss')
        axs[0].plot(self.t, vPrad/1e6, 'm--', label='Radiation+ionization')
        axs[0].plot(self.t, vPequi/1e6, 'b--', label='Equilibration')
        axs[0].plot(self.t, vPconve/1e6, 'k--', label='Electron transport')
        axs[0].plot(self.t, vPoh/1e6, 'r--', label='Ohmic heating')
        #axs[0].plot(self.t, vPnet/1e6, 'k-', label='Net electron heating power')
        axs[0].set_title('Electron energy balance')
        axs[0].legend(frameon=False)

        axs[0].set_xlim([0, self.t[-1]])
        axs[0].set_ylim([0, 0.5])

        axs[1].plot(self.t, vPequi/1e6, 'b--', label='Equilibration')
        axs[1].plot(self.t, vPcx/1e6, 'r--', label='Charge exchange')
        axs[1].plot(self.t, vPconvi/1e6, 'k--', label='Transport')
        axs[1].set_title('Ion energy balance')
        axs[1].legend(frameon=False)

        axs[1].set_xlim([0, self.t[-1]])
        axs[1].set_ylim([0, 0.5])

        fig.tight_layout()
        plt.show()


    def plotIonDensity(self, ion='D'):
        """
        Plot charge state density evolution for the named ion species.
        """
        fig, axs = plt.subplots(1,1)

        uqh = self.simulation.unknowns
        _ion = self.simulation.ions[ion]

        ni = self.x[f'ni{ion}']

        for Z0 in range(0, _ion['Z']+1):
            lbl = ion + f'{Z0}'
            if Z0 > 0: lbl += '+'

            axs.plot(self.t, ni[Z0], label=lbl)

        axs.legend(frameon=False)

        return axs


    def plotQuantity(self, data, ax=None, title=None, axis=None, *args, **kwargs):
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

        if axis is not None:
            ax.set_xlim([axis[0],axis[1]])
            ax.set_ylim([axis[2], axis[3]])

        return ax


