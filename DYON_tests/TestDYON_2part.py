# Bestäm

#'''
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

sys.path.append('../extern/DREAM/py/')
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ColdElectrons as ColdElectrons
import DREAM.Settings.Equations.ColdElectronTemperature as ColdElectronTemperature
import DREAM.Settings.Equations.ElectricField as ElectricField
import DREAM.Settings.Equations.DistributionFunction as DistFunc

sys.path.append('../py/')
from STREAM.STREAMSettings import STREAMSettings
from STREAM.STREAMOutput import STREAMOutput
from STREAM import runiface
import STREAM.Settings.Equations.IonSpecies as Ions
import STREAM.Settings.TransportSettings as Transport
#'''

# Grid parameters
#Ändra?
import numpy as np

#pMax = 1    # maximum momentum in units of m_e*c
#Np   = 300  # number of momentum grid points
#Nxi  = 20   # number of pitch grid points
tMax_initial = 5e-5  # simulation time in seconds
Nt_initial   = 500   # number of time steps
tMax_final   = 1e-2  # simulation time in seconds
Nt_final     = 1000   # number of time steps

pgp = 4.3135e-5
n_D_0 = 2.78e22 * pgp
n_D_1 = 5.56e19 * pgp
n_D = np.zeros((2,1))
#print(str(n_D))
n_D[0]=n_D_0
n_D[1]=n_D_1

n_C = 0
n_O = 0.01 * n_D_0

V_vessel = 100
B        = 2.7
a        = 0.9 # Instead of EFIT-data 0.08519 (think there's a typo should be 0.8519)
r_0      = 3   # Instead of EFIT-data 3.0381 (Should use 2.96?)
r_wall   = 1.5 # approx 1.5-2.0 m

kappa    = 1
c1       = 1.1
c2       = 0.09
c3       = 0.1

R = 7.5e-4  # Ohm, i MK2 struktur
L = 9.1e-5  # H, i MK2 struktur

r=np.array([0])
#print(str(n_D))

t   = np.linspace(0, tMax, 100)
t_d = np.array([0 , 0.02 , 0.0325, 0.0475, 0.08, 0.1 , 0.125, 0.13, 0.15, 0.20, 0.22, 0.23, 0.25, 0.3 , 0.335, 0.35, 0.37, 0.4 , 0.45, 0.5 ])
V_d = np.array([11, 21.25, 26    , 26.25 , 24  , 16.5, 8.25 , 7.9 , 7.75, 7.5 , 7.25, 6.5 , 6.5 , 6.75, 6.75 , 6   , 4.75, 4.25, 4.5 , 3.60])
V_s = interp1d(t_d, V_d, kind='linear')
V_loop_wall = V_s(t)

E_initial = V_d[0]/(2*np.pi*r_0) # Variera från 0 till V_d[0]/(2*np.pi*r_0) och se om simulering är känsligt för denna
T_e_initial = 1 # eV
T_i_initial = 0.03

sts_initial = STREAMSettings()

wall_time = L/R
print('wall_time = {} s'.format(wall_time))
sts_initial.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
sts_initial.eqsys.E_field.setInitialProfile(efield=E_initial)
sts_initial.eqsys.E_field.setBoundaryCondition(ElectricField.BC_TYPE_TRANSFORMER, V_loop_wall_R0=V_loop_wall/r_0, times=t, inverse_wall_time=1/wall_time, R0=r_0)

sts_initial.eqsys.T_cold.setType(ColdElectronTemperature.TYPE_SELFCONSISTENT)
sts_initial.eqsys.T_cold.setInitialProfile(T_e_initial)

sts_initial.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=n_D, r=r, T=T_i_initial)
sts_initial.eqsys.n_i.addIon(name='C', Z=6, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=n_C, r=r, T=T_i_initial)
sts_initial.eqsys.n_i.addIon(name='O', Z=8, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=n_O, r=r, T=T_i_initial)

sts_initial.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)
sts_initial.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_DISABLED)

sts_initial.eqsys.n_i.setJET_CWrecycling()

sts_initial.radialgrid.setB0(B)
sts_initial.radialgrid.setMinorRadius(a)
sts_initial.radialgrid.setMajorRadius(r_0)
sts_initial.radialgrid.setWallRadius(r_wall)
sts_initial.radialgrid.setVesselVolume(V_vessel)
#sts_initial.radialgrid.setElongation(kappa)
sts_initial.radialgrid.setRecyclingCoefficient1(c1)
sts_initial.radialgrid.setRecyclingCoefficient2(c2)
sts_initial.radialgrid.setRecyclingCoefficient3(c3)

sts_initial.solver.setType(Solver.NONLINEAR)
#sts_initial.solver.setDebug(savejacobian=True, savenumericaljacobian=True, timestep=1, iteration=2)
#sts_initial.solver.setDebug(savesystem=True)
sts_initial.solver.setMaxIterations(500)

#sts_initial.hottailgrid.setNxi(Nxi)
#sts_initial.hottailgrid.setNp(Np)
#sts_initial.hottailgrid.setPmax(pMax)
sts_initial.timestep.setTmax(tMax_initial)
sts_initial.timestep.setNt(Nt_initial)

sts_initial.hottailgrid.setEnabled(False)

sts_initial.runawaygrid.setEnabled(False)

sts_initial.solver.preconditioner.setEnabled(False)

sts_initial.solver.setVerbose(True)

sts_initial.other.include('fluid')

sts_initial.save('STREAMSettings_initial.h5')

sto = runiface(sts_initial, 'output_initial.h5', quiet=False)

sts_final = DREAMSettings(sts_initial)
sts_final.timestep.setTmax(tMax_final)
sts_final.timestep.setNt(Nt_final)

sts_final.fromOutput('output_initial.h5', ignore=['n_i'])
sts_final.output.setFilename('output_final.h5')
sts_final.save('STREAMSettings_final.h5')
#'''
