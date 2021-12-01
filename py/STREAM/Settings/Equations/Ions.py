# Wrapper for handling ion species

import numpy as np
import DREAM.Settings.Equations.Ions as DREAMIons
from . IonSpecies import IonSpecies, IONS_PRESCRIBED, ION_OPACITY_MODE_TRANSPARENT

from DREAM.Settings.Equations.IonSpecies import ION_CHARGED_DIFFUSION_MODE_NONE, ION_CHARGED_DIFFUSION_MODE_PRESCRIBED, ION_NEUTRAL_DIFFUSION_MODE_NONE, ION_NEUTRAL_DIFFUSION_MODE_PRESCRIBED, ION_CHARGED_ADVECTION_MODE_NONE, ION_CHARGED_ADVECTION_MODE_PRESCRIBED, ION_NEUTRAL_ADVECTION_MODE_NONE, ION_NEUTRAL_ADVECTION_MODE_PRESCRIBED


# Model to use for ion heat
IONS_T_I_NEGLECT = DREAMIons.IONS_T_I_NEGLECT
IONS_T_I_INCLUDE = DREAMIons.IONS_T_I_INCLUDE


class Ions(DREAMIons.Ions):
    

    def __init__(self, settings, *args, **kwargs):
        """
        Constructor.
        """
        super().__init__(settings=settings, *args, **kwargs)


    def addIon(self, name, Z, iontype=IONS_PRESCRIBED, Z0=None, isotope=0, SPIMolarFraction=-1, opacity_mode=ION_OPACITY_MODE_TRANSPARENT,
        charged_diffusion_mode=ION_CHARGED_DIFFUSION_MODE_NONE, charged_prescribed_diffusion=None, rChargedPrescribedDiffusion=None, tChargedPrescribedDiffusion=None,
        neutral_diffusion_mode=ION_NEUTRAL_DIFFUSION_MODE_NONE, neutral_prescribed_diffusion=None, rNeutralPrescribedDiffusion=None, tNeutralPrescribedDiffusion=None,
        charged_advection_mode=ION_CHARGED_ADVECTION_MODE_NONE, charged_prescribed_advection=None, rChargedPrescribedAdvection=None, tChargedPrescribedAdvection=None,
        neutral_advection_mode=ION_NEUTRAL_ADVECTION_MODE_NONE, neutral_prescribed_advection=None, rNeutralPrescribedAdvection=None, tNeutralPrescribedAdvection=None,
        T=None, n=None, r=None, t=None, tritium=False):
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
                         T=T, n=n, r=r, t=t, interpr=self.r, interpt=None, tritium=tritium)

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


    def setJET_CWrecycling(self, deuterium='D', carbon='C', oxygen='O', tritium='T', iron = 'Fe'):
        """
        Set recycling coefficients to the values estimated for JET-CW
        (as given in [Kim et al 2012 Nucl. Fusion 52 103016]).
        """
        for ion in self.ions:
            if ion.name == deuterium:
                ion.setRecyclingCoefficient(carbon, 0.015) #Ändrat från 0.03!!!
                ion.setRecyclingCoefficient(iron, 3e-3) # Vad ska detta vara?
            elif ion.name == oxygen:
                ion.setRecyclingCoefficient(carbon, 1.0)
                ion.setRecyclingCoefficient(oxygen, 1.0)
            elif ion.name == tritium:
                ion.setRecyclingCoefficient(tritium, 1.0)
            elif ion.name == iron:
                ion.setRecyclingCoefficient(iron, 1.0)


    def fromdict(self, data):
        """
        Load ion data from the given dictionary.
        """
        super().fromdict(data)

        # Load recycling coefficients
        if 'recycling' in data:
            rec = data['recycling'][:]
            for i in range(len(self.ions)):
                for j in range(len(self.ions)):
                    if rec[i,j] != 0.0:
                        self.ions[i].setRecyclingCoefficient(self.ions[j].name, rec[i,j])


    def todict(self):
        """
        Convert this object into a Python dictionary.
        """
        data = super().todict()

        # Construct recycling coefficient table
        nions = len(self.ions)
        # recycling coefficient table
        rec = np.zeros((nions, nions))

        for i in range(nions):
            ion = self.ions[i]

            for j in range(nions):
                rec[i,j] = ion.getRecyclingCoefficient(self.ions[j].name)

        data['recycling'] = rec

        return data


