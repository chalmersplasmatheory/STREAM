# Bestäm

#'''
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.font_manager import FontProperties

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
#import STREAM.Settings.Equations.ElectricField as ElectricField
#'''

# Grid parameters
#Ändra?
import numpy as np

#pMax = 1    # maximum momentum in units of m_e*c
#Np   = 300  # number of momentum grid points
#Nxi  = 20   # number of pitch grid points
tMax_initial = 1e-4  # simulation time in seconds
Nt_initial   = 4000   # number of time steps
tMax_final   = 5e-1  # simulation time in seconds
Nt_final     = 1e4*tMax_final   # number of time steps

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
r_wall   = 1.05 # approx 1.5-2.0 m

kappa    = 1
c1       = 1.1
c2       = 0.09
c3       = 0.1

R = 7.5e-4  # Ohm, i MK2 struktur
L = 9.1e-6  # H, i MK2 struktur

r=np.array([0])
#print(str(n_D))

#t_1   = np.linspace(0, tMax_initial, 100)
#t_2   = np.linspace(tMax_initial, 0.5, 100)
t   = np.linspace(0, 0.5, 100)
t_d = np.array([0 , 0.02 , 0.0325, 0.0475, 0.08, 0.1 , 0.125, 0.13, 0.15, 0.20, 0.22, 0.23, 0.25, 0.3 , 0.335, 0.35, 0.37, 0.4 , 0.45, 0.5 ])
V_d = np.array([11, 21.25, 26    , 26.25 , 24  , 16.5, 8.25 , 7.9 , 7.75, 7.5 , 7.25, 6.5 , 6.5 , 6.75, 6.75 , 6   , 4.75, 4.25, 4.5 , 3.60])
V_s = interp1d(t_d, V_d, kind='linear')
#V_loop_wall_1 = V_s(t_1)
#V_loop_wall_2 = V_s(t_2)
V_loop_wall = V_s(t)

#plt.plot(t,V_loop_wall)
#plt.show()

E_initial = V_loop_wall[0]/(2*np.pi*r_0) # Variera från 0 till V_d[0]/(2*np.pi*r_0) och se om simulering är känsligt för denna
E = V_loop_wall/(2*np.pi*r_0)
T_e_initial = 1 # eV
T_i_initial = 0.03
t_e=np.array([0, 0.01, 0.02,  0.03,  0.05,   0.1,   0.15,   0.2,   0.25,   0.3,   0.35,   0.4,   0.45,   0.5])
T_e=np.array([1, 2   , 7   , 10   , 42   , 152  , 206   , 250  , 277   , 294  , 312   , 320  , 330   , 335])

sts_initial = STREAMSettings()

wall_time = L/R
print('wall_time = {} s'.format(wall_time))

print(str(r_wall))
sts_initial.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
sts_initial.eqsys.E_field.setInitialProfile(efield=E_initial)
sts_initial.eqsys.E_field.setBoundaryCondition(ElectricField.BC_TYPE_TRANSFORMER, V_loop_wall_R0=V_loop_wall[0]/r_0, times=0, inverse_wall_time=1/wall_time, R0=r_0)
#sts_initial.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
#sts_initial.eqsys.E_field.setInitialProfile(efield=E_initial)
#sts_initial.eqsys.E_field.setInductances(Lp=6.09e-6, Lwall=9.1e-6, M=2.49e-6, Rwall=7.5e-4)
#sts_initial.eqsys.E_field.setCircuitVloop(V_loop_wall, t)
#sts_initial.eqsys.E_field.setPrescribedData(E, times=t)

#sts_initial.eqsys.T_cold.setType(ColdElectronTemperature.TYPE_SELFCONSISTENT)
#sts_initial.eqsys.T_cold.setInitialProfile(T_e_initial)
sts_initial.eqsys.T_cold.setPrescribedData(T_e, times=t_e)

sts_initial.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=n_D, r=r, T=T_i_initial)
#sts_initial.eqsys.n_i.addIon(name='C', Z=6, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=n_C, r=r, T=T_i_initial)
#sts_initial.eqsys.n_i.addIon(name='O', Z=8, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=n_O, r=r, T=T_i_initial)

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

#sts_initial.hottailgrid.setNxi(Nxi)
#sts_initial.hottailgrid.setNp(Np)
#sts_initial.hottailgrid.setPmax(pMax)
sts_initial.timestep.setTmax(tMax_initial)
sts_initial.timestep.setNt(Nt_initial)

sts_initial.hottailgrid.setEnabled(False)

sts_initial.runawaygrid.setEnabled(False)

sts_initial.solver.preconditioner.setEnabled(False)

#sts_initial.solver.setVerbose(True)

sts_initial.other.include('fluid', 'stream')

#sts_initial.save('STREAMSettings_initial.h5')



#sts_initial.radialgrid.setWallRadius(r_wall)

sto_initial = runiface(sts_initial, 'output_initial.h5',
                       quiet=False)
sts_final = STREAMSettings(sts_initial)
#sts_final.eqsys.E_field.setBoundaryCondition(ElectricField.BC_TYPE_TRANSFORMER, V_loop_wall_R0=V_loop_wall/r_0, times=t, inverse_wall_time=1/wall_time, R0=r_0)
sts_final.timestep.setTmax(tMax_final)
sts_final.timestep.setNt(Nt_final)
#sts_final.solver.setDebug(saveresidual=True, timestep=0, iteration=0)

sts_final.fromOutput('output_initial.h5')
sts_final.output.setFilename('output_final.h5')
sts_final.save('STREAMSettings_final.h5')

sto_final = runiface(sts_final, 'output_final.h5',
                        quiet=False)
#'''
'''
plt.plot(sto_final.grid.t[1:],sto_final.other.stream.tau_D[:],'k')
plt.plot(sto_final.grid.t[1:],sto_final.other.stream.tau_D_par[:],'r--')
plt.plot(sto_final.grid.t[1:],sto_final.other.stream.tau_D_perp[:],'b-.')
plt.legend(['tau_D','tau_D_par','tau_D_perp'])
plt.xlim(0,0.5)
plt.ylim(0,0.25)
plt.xlabel('Time [s]')
plt.ylabel('Confinement time [s]')
plt.show()

sto_final.eqsys.I_p.plot()
plt.xlim(0,0.5)
#plt.ylim(0,12e5)
plt.xlabel('Time [s]')
plt.ylabel('Plasma current [A]')
plt.show()
sto_final.eqsys.T_cold.plot()
plt.xlim(0,0.5)
#plt.ylim(0,600)
plt.xlabel('Time [s]')
plt.ylabel('Electron temperature [eV]')
plt.show()
sto_final.eqsys.n_cold.plot()
plt.xlim(0,0.5)
#plt.ylim(0,5e18)
plt.xlabel('Time [s]')
plt.ylabel('Electron density [m$^{-3}$]')
plt.show()

# Ion density
legend = []
n_D_0=sto_final.eqsys.n_i['D'][0][:]
legend.append('n_D_0')
plt.plot(sto_final.grid.t[:],n_D_0,'b-')
n_D_1=sto_final.eqsys.n_i['D'][1][:]
plt.plot(sto_final.grid.t[:],n_D_1,'b--')
legend.append('n_D_1')

n_O_0=sto_final.eqsys.n_i['O'][0][:]
legend.append('n_O_0')
plt.plot(sto_final.grid.t[:],n_O_0,'r-')
n_O_1=sto_final.eqsys.n_i['O'][1][:]
legend.append('n_O_1')
plt.plot(sto_final.grid.t[:],n_O_1,'r--')
for i in range(2,9):
    n_O_i=sto_final.eqsys.n_i['O'][i][:]
    legend.append('n_O_'+str(i))
    plt.plot(sto_final.grid.t[:],n_O_i,'r:')

n_C_0=sto_final.eqsys.n_i['C'][0][:]
legend.append('n_C_0')
plt.plot(sto_final.grid.t[:],n_C_0,'y-')
n_C_1=sto_final.eqsys.n_i['C'][1][:]
legend.append('n_C_1')
plt.plot(sto_final.grid.t[:],n_C_1,'y--')
for i in range(2,7):
    n_C_i=sto_final.eqsys.n_i['C'][i][:]
    legend.append('n_C_'+str(i))
    plt.plot(sto_final.grid.t[:],n_C_i,'y:')

plt.xlabel('Time [s]')
plt.ylabel('Ion density [m$^{-3}$]')
fontP = FontProperties()
fontP.set_size('x-small')
plt.legend(legend,prop=fontP)
#plt.xlim(0,0.12)
plt.show()
'''





# Total amount of deuterium ions
legend = []
N_D_0=sto_final.eqsys.n_i['D'][0][1:]*np.transpose(np.array(sto_final.other.stream.V_n_tot['D'][:])[np.newaxis])
legend.append('N_D_0')
plt.plot(sto_final.grid.t[1:],N_D_0,'b-')
N_D_1=sto_final.eqsys.n_i['D'][1][1:]*sto_final.other.stream.V_p[:]
plt.plot(sto_final.grid.t[1:],N_D_1,'b--')
legend.append('N_D_1')
'''
N_O_0=sto_final.eqsys.n_i['O'][0][1:]*sto_final.other.stream.V_n_tot['O'][:]
legend.append('N_O_0')
plt.plot(sto_final.grid.t[1:],N_O_0,'r-')
N_O_1=sto_final.eqsys.n_i['O'][1][1:]*sto_final.other.stream.V_p[:]
legend.append('N_O_1')
plt.plot(sto_final.grid.t[1:],N_O_1,'r--')
for i in range(2,9):
    N_O_i=sto_final.eqsys.n_i['O'][i][1:]*sto_final.other.stream.V_p[:]
    legend.append('N_O_'+str(i))
    plt.plot(sto_final.grid.t[1:],N_O_i,'r:')

N_C_0=sto_final.eqsys.n_i['C'][0][1:]*sto_final.other.stream.V_n_tot['C'][:]
legend.append('N_C_0')
plt.plot(sto_final.grid.t[1:],N_C_0,'y-')
N_C_1=sto_final.eqsys.n_i['C'][1][1:]*sto_final.other.stream.V_p[:]
legend.append('N_C_1')
plt.plot(sto_final.grid.t[1:],N_C_1,'y--')
for i in range(2,7):
    N_C_i=sto_final.eqsys.n_i['C'][i][1:]*sto_final.other.stream.V_p[:]
    legend.append('N_C_'+str(i))
    plt.plot(sto_final.grid.t[1:],N_C_i,'y:')
'''
plt.xlabel('Time [s]')
plt.ylabel('Total amount of ions per state')
fontP = FontProperties()
fontP.set_size('x-small')
plt.legend(legend,prop=fontP)
#plt.xlim(0,0.12)
plt.show()

# Total amount of deuterium ions
legend = []
N_D_0=sto_final.eqsys.n_i['D'][0][1:]*np.transpose(np.array(sto_final.other.stream.V_n_tot['D'][:])[np.newaxis])
N_D_1=sto_final.eqsys.n_i['D'][1][1:]*sto_final.other.stream.V_p[:]
N_D=N_D_0+N_D_1
plt.plot(sto_final.grid.t[1:],N_D,'b-')
legend.append('N_D')
'''
N_O_0=sto_final.eqsys.n_i['O'][0][1:]*sto_final.other.stream.V_n_tot['O'][:]
N_O_1=sto_final.eqsys.n_i['O'][1][1:]*sto_final.other.stream.V_p[:]
N_O=N_O_0+N_O_1
for i in range(2,9):
    N_O+=sto_final.eqsys.n_i['O'][i][1:]*sto_final.other.stream.V_p[:]
legend.append('N_O')
plt.plot(sto_final.grid.t[1:],N_O,'r-')

N_C_0=sto_final.eqsys.n_i['C'][0][1:]*sto_final.other.stream.V_n_tot['C'][:]
N_C_1=sto_final.eqsys.n_i['C'][1][1:]*sto_final.other.stream.V_p[:]
N_C = N_C_0 + N_C_1
for i in range(2,7):
    N_C+=sto_final.eqsys.n_i['C'][i][1:]*sto_final.other.stream.V_p[:]
legend.append('N_C')
plt.plot(sto_final.grid.t[1:],N_C,'y-')
'''
plt.xlabel('Time [s]')
plt.ylabel('Total amount of ions per species')
fontP = FontProperties()
fontP.set_size('x-small')
plt.legend(legend,prop=fontP)
#plt.xlim(0,0.12)
plt.show()

# Total amount of ions
N_D_0=sto_final.eqsys.n_i['D'][0][1:]*np.transpose(np.array(sto_final.other.stream.V_n_tot['D'][:])[np.newaxis])
N_D_1=sto_final.eqsys.n_i['D'][1][1:]*sto_final.other.stream.V_p[:]
N= N_D_0+N_D_1
'''
N_O_0=sto_final.eqsys.n_i['O'][0][1:]*sto_final.other.stream.V_n_tot['O'][:]
N_O_1=sto_final.eqsys.n_i['O'][1][1:]*sto_final.other.stream.V_p[:]
N+=N_O_0+N_O_1
for i in range(2,9):
    N+=sto_final.eqsys.n_i['O'][i][1:]*sto_final.other.stream.V_p[:]

N_C_0=sto_final.eqsys.n_i['C'][0][1:]*sto_final.other.stream.V_n_tot['C'][:]
N_C_1=sto_final.eqsys.n_i['C'][1][1:]*sto_final.other.stream.V_p[:]
N += N_C_0 + N_C_1
for i in range(2,7):
    N+=sto_final.eqsys.n_i['C'][i][1:]*sto_final.other.stream.V_p[:]
'''
plt.plot(sto_final.grid.t[1:],N)

plt.xlabel('Time [s]')
plt.ylabel('Total amount of ions')
fontP = FontProperties()
fontP.set_size('x-small')
#plt.legend(legend,prop=fontP)
#plt.xlim(0,0.12)
plt.ylim(bottom=0)
plt.show()
#'''
