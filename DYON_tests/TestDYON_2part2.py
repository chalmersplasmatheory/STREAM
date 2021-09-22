# Bestäm

#'''
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.constants
from matplotlib.font_manager import FontProperties
from PlasmaParameters import evaluateSpitzerConductivity

sys.path.append('../extern/DREAM/py/')
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ColdElectrons as ColdElectrons
import DREAM.Settings.Equations.ColdElectronTemperature as ColdElectronTemperature
import DREAM.Settings.Equations.ElectricField as ElectricField_D
import DREAM.Settings.Equations.DistributionFunction as DistFunc

sys.path.append('../py/')
from STREAM.STREAMSettings import STREAMSettings
from STREAM.STREAMOutput import STREAMOutput
from STREAM import runiface
import STREAM.Settings.Equations.IonSpecies as Ions
import STREAM.Settings.TransportSettings as Transport
import STREAM.Settings.Equations.ElectricField as ElectricField
#'''

# Grid parameters
#Ändra?

#pMax = 1    # maximum momentum in units of m_e*c
#Np   = 300  # number of momentum grid points
#Nxi  = 20   # number of pitch grid points
tMax_initial = 1e-4  # simulation time in seconds
Nt_initial   = 2000   # number of time steps
tMax_final   = 3e-1-tMax_initial  # simulation time in seconds
Nt_final     = 15e3*tMax_final   # number of time steps

pgp = 2.7e-3
gamma_i = 2e-3      # Ionization fraction
n_D_tot = 4.8e20 * pgp
n_D = np.zeros((2,1))
#print(str(n_D))
n_D[0]=n_D_tot*(1-gamma_i)
n_D[1]=n_D_tot*gamma_i

n_C = 0
n_O = 0.001 * n_D_tot

V_vessel = 100
B        = 2.4
r_0      = 2.96
r_wall   = 1

t   = np.linspace(0, tMax_final, 100)
t_a   = np.array([0  ,  0.017,  0.05,  0.085,  0.14,  0.19,  0.25,  0.3])
V_p   = np.array([100, 80    , 56   , 48    , 52   , 51.75, 54.25, 56  ])
a_vec = np.sqrt(V_p/(2*np.pi**2*r_0))
a_fun = interp1d(t_a, a_vec, kind='cubic')
a = a_fun(t)


kappa    = 1
c1       = 1.1
c2       = 0.05
c3       = 0.1

R = 7.5e-4  # Ohm, i MK2 struktur
L = 9.1e-6  # H, i MK2 struktur
Lp = scipy.constants.mu_0*r_0*(np.log(8*r_0/a[0]) + 0.5 - 2)
print(str(Lp))

r=np.array([0])

t_Vloop = [0 , 0.02 , 0.0325, 0.0475, 0.08, 0.1 , 0.125, 0.13, 0.15, 0.20, 0.22, 0.23, 0.25, 0.3 , 0.335, 0.35, 0.37, 0.4 , 0.45, 0.5 ]
d_Vloop = [11, 21.25, 26    , 26.25 , 24  , 16.5, 8.25 , 7.9 , 7.75, 7.5 , 7.25, 6.5 , 6.5 , 6.75, 6.75 , 6   , 4.75, 4.25, 4.5 , 3.60]
V_s = interp1d(t_Vloop, d_Vloop, kind='linear')
V_loop_wall = V_s(t)

T_e_initial = 1
T_i_initial = 0.026
t_e=np.array([0, 0.01, 0.02,  0.03,  0.05,   0.1,   0.15,   0.2,   0.25,   0.3,   0.35,   0.4,   0.45,   0.5])
T_e=np.array([1, 2   , 7   , 10   , 42   , 152  , 206   , 250  , 277   , 294  , 312   , 320  , 330   , 335])

sigma = evaluateSpitzerConductivity(n=n_D[1], T=T_e_initial, Z=1)
I_p= 2.4e3
J = I_p/(a_vec[0]**2 *np.pi) # 943?
E_initial = J/sigma

sts_initial = STREAMSettings()

wall_time = L/R*0.3
print('wall_time = {} s'.format(wall_time))

print(str(r_wall))
#sts_initial.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
#sts_initial.eqsys.E_field.setInitialProfile(efield=E_initial)
#sts_initial.eqsys.E_field.setBoundaryCondition(ElectricField_D.BC_TYPE_TRANSFORMER, V_loop_wall_R0=V_loop_wall/r_0, times=t, inverse_wall_time=1/wall_time, R0=r_0)
sts_initial.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
sts_initial.eqsys.E_field.setInitialProfile(efield=E_initial)
sts_initial.eqsys.E_field.setInductances(Lp=5.19e-6, Lwall=9.1e-6, M=2.49e-6, Rwall=7.5e-4)
sts_initial.eqsys.E_field.setCircuitVloop(V_loop_wall, t)

sts_initial.eqsys.T_cold.setType(ColdElectronTemperature.TYPE_SELFCONSISTENT)
sts_initial.eqsys.T_cold.setInitialProfile(T_e_initial)
#sts_initial.eqsys.T_cold.setPrescribedData(T_e, times=t_e)

sts_initial.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=n_D, r=r, T=T_i_initial)
sts_initial.eqsys.n_i.addIon(name='C', Z=6, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=n_C, r=r, T=T_i_initial)
sts_initial.eqsys.n_i.addIon(name='O', Z=8, iontype=Ions.IONS_DYNAMIC_NEUTRAL, n=n_O, r=r, T=T_i_initial)

sts_initial.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)
sts_initial.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_DISABLED)

sts_initial.eqsys.n_i.setJET_CWrecycling()

sts_initial.radialgrid.setB0(B)
sts_initial.radialgrid.setMinorRadius(a, t=t)
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

sts_initial.other.include('fluid', 'stream', 'scalar')

sts_initial.save('settings_initial.h5')



#sts_initial.radialgrid.setWallRadius(r_wall)

sto_initial = runiface(sts_initial, 'output_initial.h5',
                       quiet=False)
sts_final = STREAMSettings(sts_initial)
#sts_final.eqsys.E_field.setBoundaryCondition(ElectricField.BC_TYPE_TRANSFORMER, V_loop_wall_R0=V_loop_wall/r_0, times=t, inverse_wall_time=1/wall_time, R0=r_0)
sts_final.timestep.setTmax(tMax_final)
sts_final.timestep.setNt(Nt_final)
sts_final.timestep.setNumberOfSaveSteps(1000)
#sts_final.solver.setDebug(saveresidual=True, timestep=0, iteration=0)

sts_final.radialgrid.setMinorRadius(a, t=t-sto_initial.grid.t[-1])

sts_final.fromOutput('output_initial.h5')
sts_final.output.setFilename('output_final.h5')
sts_final.save('settings_final.h5')

sto_final = runiface(sts_final, 'output_final.h5',
                        quiet=False)
#'''

#print(str(sto_initial.eqsys.I_p[:]))
fig, axs = plt.subplots(3,2)

# Plasma current
axs[0,0].plot(sto_final.grid.t[:],sto_final.eqsys.I_p[:])
axs[0,0].set_xlim([0,0.3])
axs[0,0].set_ylim([0,6e5])
axs[0,0].set_xlabel('Time [s]')
axs[0,0].set_ylabel('Plasma current [A]')

# Electron temperature
axs[1,0].plot(sto_final.grid.t[:],sto_final.eqsys.T_cold[:])
axs[1,0].set_xlim([0,0.3])
axs[1,0].set_ylim([0,400])
axs[1,0].set_xlabel('Time [s]')
axs[1,0].set_ylabel('Electron temperature [eV]')

# Ion temperature
axs[2,0].plot(sto_final.grid.t[:],2.0/3.0*sto_final.eqsys.W_i['D'][:]/sto_final.eqsys.N_i['D'][:]/1.60217662e-19)
axs[2,0].set_xlim([0,0.3])
axs[2,0].set_ylim([0,400])
axs[2,0].set_xlabel('Time [s]')
axs[2,0].set_ylabel('Ion temperature [eV]')

# Electron density
axs[0,1].plot(sto_final.grid.t[:],sto_final.eqsys.n_cold[:])
axs[0,1].set_xlim([0,0.3])
axs[0,1].set_ylim([0,6e18])
axs[0,1].set_xlabel('Time [s]')
axs[0,1].set_ylabel('Electron density [m$^{-3}$]')

# Confinement time
axs[2,1].plot(sto_final.grid.t[1:],sto_final.other.stream.tau_D[:])
axs[2,1].set_xlim([0,0.3])
axs[2,1].set_ylim([0,0.1])
axs[2,1].set_xlabel('Time [s]')
axs[2,1].set_ylabel('Confinement time [s]')

plt.show()
'''
# Ion density
legend = []
n_D_0=sto_final.eqsys.n_i['D'][0][:]
legend.append('n_D_0')
plt.plot(sto_final.grid.t[:],n_D_0,'b-')
n_D_1=sto_final.eqsys.n_i['D'][1][:]
plt.plot(sto_final.grid.t[:],n_D_1,'b--')
legend.append('n_D_1')
#''' '''
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
#''' '''
plt.xlabel('Time [s]')
plt.ylabel('Ion density [m$^{-3}$]')
fontP = FontProperties()
fontP.set_size('x-small')
plt.legend(legend,prop=fontP)
#plt.xlim(0,0.12)
plt.show()

n_D_0 = np.array(n_D_0).flatten()
dnD0_idt = np.gradient(n_D_0)
#''' '''
n_O_0 = np.array(n_O_0).flatten()
dnO0_idt = np.gradient(n_O_0)
n_C_0 = np.array(n_C_0).flatten()
dnC0_idt = np.gradient(n_C_0)
#''' '''
plt.plot(sto_final.grid.t[:], dnD0_idt)
plt.plot(sto_final.grid.t[:], dnO0_idt)
plt.plot(sto_final.grid.t[:], dnC0_idt)
plt.plot(np.linspace(0,tMax_final,2),np.linspace(0,0,2),'k')
plt.xlabel('Time [s]')
plt.ylabel('Derivative of ion density [m$^{-3}$]')
#fontP = FontProperties()
#fontP.set_size('x-small')
plt.legend(['dn_D_0/dt','dn_O_0/dt','dn_C_0/dt'],prop=fontP)
#plt.xlim(0,0.12)
plt.show()


# Total amount of deuterium ions
legend = []
N_D_0=sto_final.eqsys.n_i['D'][0][1:]*sto_final.other.stream.V_n_tot['D'][:]#np.transpose(np.array([np.newaxis]))
legend.append('N_D_0')
plt.plot(sto_final.grid.t[1:],N_D_0,'b-')
N_D_1=sto_final.eqsys.n_i['D'][1][1:]*sto_final.other.stream.V_p[:]
plt.plot(sto_final.grid.t[1:],N_D_1,'b--')
legend.append('N_D_1')
#''' '''
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
#''' '''
plt.xlabel('Time [s]')
plt.ylabel('Total amount of ions per state')
fontP = FontProperties()
fontP.set_size('x-small')
plt.legend(legend,prop=fontP)
#plt.xlim(0,0.12)
plt.show()

# Total amount of deuterium ions
legend = []
N_D_0=sto_final.eqsys.n_i['D'][0][1:]*sto_final.other.stream.V_n_tot['D'][:]#np.transpose(np.array([np.newaxis]))
N_D_1=sto_final.eqsys.n_i['D'][1][1:]*sto_final.other.stream.V_p[:]
N_D=N_D_0+N_D_1
plt.plot(sto_final.grid.t[1:],N_D,'b-')
legend.append('N_D')
#''' '''
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
#''' '''
plt.xlabel('Time [s]')
plt.ylabel('Total amount of ions per species')
fontP = FontProperties()
fontP.set_size('x-small')
plt.legend(legend,prop=fontP)
#plt.xlim(0,0.12)
plt.show()

# Total amount of ions
N_D_0=sto_final.eqsys.n_i['D'][0][1:]*sto_final.other.stream.V_n_tot['D'][:]#np.transpose(np.array([np.newaxis]))
N_D_1=sto_final.eqsys.n_i['D'][1][1:]*sto_final.other.stream.V_p[:]
N= N_D_0+N_D_1
#''' '''
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
#''' '''
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
