# Wrapper for handling ion species

import numpy as np
import DREAM.Settings.Equations.Ions as DREAMIons
from . IonSpecies import IonSpecies, IONS_PRESCRIBED, ION_OPACITY_MODE_TRANSPARENT

from DREAM.Settings.Equations.IonSpecies import ION_CHARGED_DIFFUSION_MODE_NONE, ION_CHARGED_DIFFUSION_MODE_PRESCRIBED, ION_NEUTRAL_DIFFUSION_MODE_NONE, ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED, ION_CHARGED_ADVECTION_MODE_NONE, ION_CHARGED_ADVECTION_MODE_PRESCRIBED, ION_NEUTRAL_ADVECTION_MODE_NONE, ION_NEUTRAL_ADVECTION_MODE_PRESCRIBED
from DREAM.Settings.Equations.PrescribedParameter import PrescribedParameter


# Model to use for ion heat
IONS_T_I_NEGLECT = DREAMIons.IONS_T_I_NEGLECT
IONS_T_I_INCLUDE = DREAMIons.IONS_T_I_INCLUDE


class Ions(DREAMIons.Ions, PrescribedParameter):
    

    def __init__(self, settings, *args, **kwargs):
        """
        Constructor.
        """
        super().__init__(settings=settings, *args, **kwargs)

        self.fueling = {}


    def addIon(self, name, Z, iontype=IONS_PRESCRIBED, Z0=None, isotope=0, SPIMolarFraction=-1, opacity_mode=ION_OPACITY_MODE_TRANSPARENT,
        charged_diffusion_mode=ION_CHARGED_DIFFUSION_MODE_NONE, charged_prescribed_diffusion=None, rChargedPrescribedDiffusion=None, tChargedPrescribedDiffusion=None,
        neutral_diffusion_mode=ION_NEUTRAL_DIFFUSION_MODE_NONE, neutral_prescribed_diffusion=None, rNeutralPrescribedDiffusion=None, tNeutralPrescribedDiffusion=None,
        charged_advection_mode=ION_CHARGED_ADVECTION_MODE_NONE, charged_prescribed_advection=None, rChargedPrescribedAdvection=None, tChargedPrescribedAdvection=None,
        neutral_advection_mode=ION_NEUTRAL_ADVECTION_MODE_NONE, neutral_prescribed_advection=None, rNeutralPrescribedAdvection=None, tNeutralPrescribedAdvection=None,
        T=None, n=None, r=None, t=None, tritium=False, hydrogen=False):
        """
        Adds a new ion species to the plasma.

        :param str name:               Name by which the ion species will be referred to.
        :param int Z:                  Ion charge number.
        :param int isotope:            Ion mass number.
        :param int iontype:            Method to use for evolving ions in time.
        :param int Z0:                 Charge state to populate (used for populating exactly one charge state for the ion).
        :param n:                      Ion density (can be either a scalar, 1D array or 2D array, depending on the other input parameters)
        :param float SPIMolarFraction: Molar fraction of the SPI injection (if any). A negative value means that this species is not part of the SPI injection 
        :param numpy.ndarray r:        Radial grid on which the input density is defined.
        :param T:                      Ion initial temperature (can be scalar for uniform temperature, otherwise 1D array matching `r` in size)
        :param numpy.ndarray r:        Radial grid on which the input density and temperature is defined.
        :param numpy.ndarray t:        Time grid on which the input density is defined.
        :param bool tritium:           If ``True``, the ion species is treated as Tritium.
        """
        if (self.r is not None) and (r is not None) and (np.any(self.r != r)):
            raise EquationException("The radial grid must be the same for all ion species.")
        if (self.t is not None) and (t is not None) and (np.any(self.t != t)):
            raise EquationException("The time grid must be the same for all ion species.")
        if (self.rChargedPrescribedDiffusion is not None) and (rChargedPrescribedDiffusion is not None) and (
        np.any(self.rChargedPrescribedDiffusion != rChargedPrescribedDiffusion)):
            raise EquationException(
                "The radial grid for the prescribed charged diffusion must be the same for all ion species.")
        if (self.tChargedPrescribedDiffusion is not None) and (tChargedPrescribedDiffusion is not None) and (
        np.any(self.tChargedPrescribedDiffusion != tChargedPrescribedDiffusion)):
            raise EquationException(
                "The time grid for the prescribed charged diffusion must be the same for all ion species.")

        if (self.rNeutralPrescribedDiffusion is not None) and (rNeutralPrescribedDiffusion is not None) and (
        np.any(self.rNeutralPrescribedDiffusion != rNeutralPrescribedDiffusion)):
            raise EquationException(
                "The radial grid for the prescribed neutral diffusion must be the same for all ion species.")
        if (self.tNeutralPrescribedDiffusion is not None) and (tNeutralPrescribedDiffusion is not None) and (
        np.any(self.tNeutralPrescribedDiffusion != tNeutralPrescribedDiffusion)):
            raise EquationException(
                "The time grid for the prescribed neutral diffusion must be the same for all ion species.")

        if T is not None:
            self.typeTi = IONS_T_I_INCLUDE

        ion = IonSpecies(settings=self.settings, name=name, Z=Z, ttype=iontype, Z0=Z0, isotope=isotope, SPIMolarFraction=SPIMolarFraction, opacity_mode=opacity_mode,
                         charged_diffusion_mode=charged_diffusion_mode,
                         charged_prescribed_diffusion=charged_prescribed_diffusion,
                         rChargedPrescribedDiffusion=rChargedPrescribedDiffusion,
                         tChargedPrescribedDiffusion=tChargedPrescribedDiffusion,
                         neutral_diffusion_mode=neutral_diffusion_mode,
                         neutral_prescribed_diffusion=neutral_prescribed_diffusion,
                         rNeutralPrescribedDiffusion=rNeutralPrescribedDiffusion,
                         tNeutralPrescribedDiffusion=tNeutralPrescribedDiffusion,
                         charged_advection_mode=charged_advection_mode,
                         charged_prescribed_advection=charged_prescribed_advection,
                         rChargedPrescribedAdvection=rChargedPrescribedAdvection,
                         tChargedPrescribedAdvection=tChargedPrescribedAdvection,
                         neutral_advection_mode=neutral_advection_mode,
                         neutral_prescribed_advection=neutral_prescribed_advection,
                         rNeutralPrescribedAdvection=rNeutralPrescribedAdvection,
                         tNeutralPrescribedAdvection=tNeutralPrescribedAdvection,
                         T=T, n=n, r=r, t=t, interpr=self.r, interpt=None, tritium=tritium, hydrogen=hydrogen)

        self.ions.append(ion)

        self.r = ion.getR()
        if ion.getTime() is not None:
            self.t = ion.getTime()

        if charged_diffusion_mode == ION_CHARGED_DIFFUSION_MODE_PRESCRIBED:
            self.rChargedPrescribedDiffusion = ion.getRChargedPrescribedDiffusion()
            self.tChargedPrescribedDiffusion = ion.getTChargedPrescribedDiffusion()

        if neutral_diffusion_mode == ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED:
            self.rNeutralPrescribedDiffusion = ion.getRNeutralPrescribedDiffusion()
            self.tNeutralPrescribedDiffusion = ion.getTNeutralPrescribedDiffusion()

        if charged_advection_mode == ION_CHARGED_ADVECTION_MODE_PRESCRIBED:
            self.rChargedPrescribedAdvection = ion.getRChargedPrescribedAdvection()
            self.tChargedPrescribedAdvection = ion.getTChargedPrescribedAdvection()

        if neutral_advection_mode == ION_NEUTRAL_ADVECTION_MODE_PRESCRIBED:
            self.rNeutralPrescribedAdvection = ion.getRNeutralPrescribedAdvection()
            self.tNeutralPrescribedAdvection = ion.getTNeutralPrescribedAdvection()

    def getIndex(self, species):
        for iIon in range(len(self.ions)):
            if self.ions[iIon].getName() == species:
                return iIon
        print("WARNING: no species of name '{}'.".format(species))
        return None


    def setJET_CWrecycling(self, deuterium='D', carbon='C', oxygen='O', tritium='T', Nt=100, tMax=20):
        """
        Set recycling coefficients to the values estimated for JET-CW
        (as given in [Kim et al 2012 Nucl. Fusion 52 103016]).
        """
        t = np.linspace(0, tMax, Nt)
        for ion in self.ions:
            if ion.name == deuterium:
                c1 = 1.1
                c2 = 0.05
                c3 = 0.1
                Y = c1 - c2 * (1 - np.exp(-t/c3))
                ion.setRecyclingCoefficient(deuterium, Y, t)
                ion.setRecyclingCoefficient(carbon, 0.015)
            elif ion.name == oxygen:
                ion.setRecyclingCoefficient(carbon, 1.0)
                ion.setRecyclingCoefficient(oxygen, 1.0)
            elif ion.name == tritium:
                c1 = 1.1
                c2 = 0.05
                c3 = 0.1
                Y = c1 - c2 * (1 - np.exp(-t / c3))
                ion.setRecyclingCoefficient(tritium, Y, t)


    def setFueling(self, species, fueling, times=0):
        """
        Prescribes the time evolution of the fueling source for the specified
        ion species.
        """
        n, r, t = self._setPrescribedData(data=fueling, times=times)

        if len(self.fueling) > 0:
            k = list(self.fueling.keys())[0]
            if r != self.fueling[k]['r']:
                raise EquationException("ions: All fueling functions must have the same 'r' grids.")
            if (t != self.fueling[k]['t']).any():
                raise EquationException("ions: All fueling functions must have the same 't' grids.")

        self.fueling[species] = {
            'n': n,
            't': t,
            'r': r
        }


    def fromdict(self, data):
        """
        Load ion data from the given dictionary.
        """
        super().fromdict(data)

        names = data['names'].split(';')[:-1]
        Z = data['Z']
        types = data['types']

        # Load recycling coefficients
        if 'recycling' in data:
            sputRecCoefficientTable = data['recycling']['x'][:]
            tSputRec = data['recycling']['t'][:]
            for i in range(len(self.ions)):
                for j in range(len(self.ions)):
                    self.ions[i].setRecyclingCoefficient(self.ions[j].name, sputRecCoefficientTable[i,j,:], tSputRec)

        if 'fueling' in data:
            idx = 0
            for i in range(len(Z)):
                if types[i] == IONS_PRESCRIBED:
                    continue

                n = data['fueling']['x'][idx,:]

                idx += Z[i]+1
                if np.all(n==0):
                    continue

                self.fueling[names[i]] = {
                    'n': n,
                    'r': data['fueling']['r'],
                    't': data['fueling']['t']
                }


    def todict(self):
        """
        Convert this object into a Python dictionary.
        """
        data = super().todict()

        # Construct recycling coefficient table
        nions = len(self.ions)

        # recycling coefficient table and corresponding time array
        breakflag = False
        for i in range(nions):
            for j in range(nions):
                tSputRec = self.ions[i].getRecyclingCoefficient(self.ions[j].name)[1]
                NtSputRec = len(tSputRec)
                if NtSputRec != 1:
                    breakflag = True
                    break
            if breakflag:
                break


        sputRecCoefficientTable = np.zeros((nions, nions, NtSputRec))
        for i in range(nions):
            ion = self.ions[i]
            for j in range(nions):
                values, tvalues = ion.getRecyclingCoefficient(self.ions[j].name)
                if len(values) == 1:
                    sputRecCoefficientTable[i,j,:] = np.array(values[0]*np.ones(NtSputRec))
                elif any(tvalues != tSputRec):
                    raise EquationException(
                        "Ions: Recycling coefficient for species '{}' due to species '{}' does not have the correct time vector.".format(
                            self.ions[j], ion))
                elif len(values) != NtSputRec:
                    raise EquationException(
                        "Ions: Recycling coefficient for species '{}' due to species '{}' does not have the correct length.".format(
                            self.ions[j], ion))
                else:
                    sputRecCoefficientTable[i,j,:] = values

        data['recycling'] = {
            'x': sputRecCoefficientTable,
            't': tSputRec
        }

        if len(self.fueling) > 0:
            # Calculate number of charge states
            nZ0 = 0
            for ion in self.ions:
                if ion.ttype == IONS_PRESCRIBED:
                    continue
                nZ0 += ion.Z+1

            k = list(self.fueling.keys())[0]
            t, r   = self.fueling[k]['t'], self.fueling[k]['r']
            nt, nr = t.size, r.size
            fueling = np.zeros((nZ0, nt, nr))
            idx = 0
            for ion in self.ions:
                if ion.ttype == IONS_PRESCRIBED:
                    continue

                if ion.name in self.fueling:
                    # Only set the neutral state...
                    fueling[idx,:] = self.fueling[ion.name]['n']


                idx += ion.Z+1
            data['fueling'] = {
                'x': fueling,
                'r': r,
                't': t
            }

        return data


