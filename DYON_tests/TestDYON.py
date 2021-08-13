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
import DREAM.Settings.Equations.DistributionFunction as DistFunc

sys.path.append('../py/')
from STREAM.STREAMSettings import STREAMSettings
from STREAM.STREAMOutput import STREAMOutput
from STREAM import runiface
import STREAM.Settings.Equations.ElectricField as ElectricField
import STREAM.Settings.Equations.IonSpecies as Ions
import STREAM.Settings.TransportSettings as Transport
#'''

# Grid parameters
#Ändra?
import numpy as np

#pMax = 1    # maximum momentum in units of m_e*c
#Np   = 300  # number of momentum grid points
#Nxi  = 20   # number of pitch grid points
tMax = 5e-1  # simulation time in seconds
Nt   = 400000 # number of time steps

pgp = 4.3135e-5
n_D_0 = 2.78e22 * pgp
n_D_1 = 5.56e19 * pgp
n_D = np.zeros((2,1))
#print(str(n_D))
n_D[0]=n_D_0
n_D[1]=n_D_1

n_C = 0
n_O = 0.01 * n_D_0

V_vessel = 50#100
B        = 2.7
a        = 0.9 # Instead of EFIT-data 0.08519 (think there's a typo should be 0.8519)
r_0      = 3   # Instead of EFIT-data 3.0381 (Should use 2.96?)
#r_0      = 1.2
r_wall   = 1.5 #1/np.pi*np.sqrt(V_vessel/(2*r_0))

kappa    = 1
c1       = 1.1
c2       = 0.09
c3       = 0.1

R = 7.5e-4  # Ohm, i MK2 struktur
L = 9.1e-5  # H, i MK2 struktur

r=np.array([0])
#print(str(n_D))

t   = np.linspace(0, 0.5, 100)
t_d = np.array([0 , 0.02 , 0.0325, 0.0475, 0.08, 0.1 , 0.125, 0.13, 0.15, 0.20, 0.22, 0.23, 0.25, 0.3 , 0.335, 0.35, 0.37, 0.4 , 0.45, 0.5 ])
V_d = np.array([11, 21.25, 26    , 26.25 , 24  , 16.5, 8.25 , 7.9 , 7.75, 7.5 , 7.25, 6.5 , 6.5 , 6.75, 6.75 , 6   , 4.75, 4.25, 4.5 , 3.60])
V_s = interp1d(t_d, V_d, kind='linear')
V_loop_wall = V_s(t)

E_initial = 0.04 #V_d[0]/(2*np.pi*r_0) # Variera från 0 till V_d[0]/(2*np.pi*r_0) och se om simulering är känsligt för denna
T_e_initial = 1 # eV
T_i_initial = 0.03

t_e=np.array([0, 0.01, 0.02,  0.03,  0.05,   0.1,   0.15 ,  0.2,   0.25,   0.3 ,  0.35,   0.4,   0.45,   0.5])
T_e=np.array([1, 2,    7 ,   10,    42,    152,  206,    250,   277,    294,   312,    320,   330,    335])

sts = STREAMSettings()

"""
wall_time = L/R
print('wall_time = {} s'.format(wall_time))
sts.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
sts.eqsys.E_field.setInitialProfile(efield=E_initial)
sts.eqsys.E_field.setBoundaryCondition(ElectricField.BC_TYPE_TRANSFORMER, V_loop_wall_R0=V_loop_wall/r_0, times=t, inverse_wall_time=1/wall_time, R0=r_0)
"""
sts.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
sts.eqsys.E_field.setInitialProfile(efield=E_initial)
sts.eqsys.E_field.setInductances(Lp=5.4e-6, Lwall=L, M=2.49e-6, Rwall=7.5e-4)
sts.eqsys.E_field.setCircuitVloop(V_loop_wall, t)

#sts.eqsys.T_cold.setType(ColdElectronTemperature.TYPE_SELFCONSISTENT)
#sts.eqsys.T_cold.setInitialProfile(T_e_initial)
sts.eqsys.T_cold.setPrescribedData(T_e, times=t_e)

sts.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=n_D, r=r, T=T_i_initial)
sts.eqsys.n_i.addIon(name='C', Z=6, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=n_C, r=r, T=T_i_initial)
sts.eqsys.n_i.addIon(name='O', Z=8, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=n_O, r=r, T=T_i_initial)

sts.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)
sts.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_DISABLED)

sts.eqsys.n_i.setJET_CWrecycling()

sts.radialgrid.setB0(B)
sts.radialgrid.setMinorRadius(a)
sts.radialgrid.setMajorRadius(r_0)
sts.radialgrid.setWallRadius(r_wall)
sts.radialgrid.setVesselVolume(V_vessel)
#sts.radialgrid.setElongation(kappa)
sts.radialgrid.setRecyclingCoefficient1(c1)
sts.radialgrid.setRecyclingCoefficient2(c2)
sts.radialgrid.setRecyclingCoefficient3(c3)

sts.solver.setType(Solver.NONLINEAR)
#sts.solver.setDebug(savejacobian=True, savenumericaljacobian=True, timestep=1, iteration=2)
#sts.solver.setDebug(savesystem=True)

#sts.hottailgrid.setNxi(Nxi)
#sts.hottailgrid.setNp(Np)
#sts.hottailgrid.setPmax(pMax)
sts.timestep.setTmax(tMax)
sts.timestep.setNt(Nt)
sts.timestep.setNumberOfSaveSteps(Nt/100)

sts.hottailgrid.setEnabled(False)

sts.runawaygrid.setEnabled(False)

sts.solver.preconditioner.setEnabled(False)
#sts.solver.tolerance.set('lambda_i', reltol=1e-4)
#sts.solver.tolerance.set('n_i', reltol=1e-4)
#sts.solver.tolerance.set('n_tot', reltol=1e-4)
#sts.solver.tolerance.set('N_i', reltol=1e-4)

#sts.solver.setVerbose(True)

#sts.solver.setDebug(printjacobianinfo=True, timestep=1, iteration=0, savesystem=True)

sts.other.include('fluid', 'stream')

sts.save('STREAMSettings.h5')

sto = runiface(sts, 'output.h5', quiet=False)

'''
sts2 = STREAMSettings(sts)
sts2.timestep.setTmax(tMax)
sts2.timestep.setNt(Nt/10)

sto2 = runiface(sts2, 'output2.h5', quiet=False)
#'''
