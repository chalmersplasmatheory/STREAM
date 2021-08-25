#'''
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.font_manager import FontProperties
from PlasmaParameters import evaluateSpitzerConductivity
import scipy.constants

sys.path.append('../extern/DREAM/py/')
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ColdElectrons as ColdElectrons
import DREAM.Settings.Equations.ColdElectronTemperature as ColdElectronTemperature
import DREAM.Settings.Equations.ColdElectronTemperature as Tcold
import DREAM.Settings.Equations.ElectricField as ElectricField
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
import numpy as np

#pMax = 1    # maximum momentum in units of m_e*c
#Np   = 300  # number of momentum grid points
#Nxi  = 20   # number of pitch grid points
tMax_initial = 1e-4  # simulation time in seconds
Nt_initial   = 4000   # number of time steps
tMax_final   = 3e-2  # simulation time in seconds
Nt_final     = 1e4*tMax_final   # number of time steps

gamma=2e-3
prefill=5e-5
n0 = 3.22e22 * prefill  # Initial total deuterium density
nD = n0 * np.array([[1-gamma], [gamma]])

V_vessel = 100
B        = 2.65
a        = 1.6
r_0      = 5.65
r_wall   = 1

Vloop=12

Lp = float(scipy.constants.mu_0 * r_0 * (np.log(8*r_0/a) + 0.5 - 2))

r=np.array([0])

T_e_initial = 1
T_i_initial = 0.026

a_T=0.5


sigma = evaluateSpitzerConductivity(n=nD[1], T=T_e_initial, Z=1)
J = 883.3 # 943?
E_initial = J/sigma

sto_list=[]

sts_initial = STREAMSettings()

sts_initial.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
sts_initial.eqsys.E_field.setInitialProfile(efield=E_initial)
sts_initial.eqsys.E_field.setInductances(Lp=Lp, Lwall=9.1e-6, M=2.49e-6, Rwall=1e6)
sts_initial.eqsys.E_field.setCircuitVloop(Vloop)

sts_initial.eqsys.T_cold.setType(ColdElectronTemperature.TYPE_SELFCONSISTENT)
sts_initial.eqsys.T_cold.setInitialProfile(T_e_initial)
sts_initial.eqsys.T_cold.setRecombinationRadiation(recombination=Tcold.RECOMBINATION_RADIATION_INCLUDED)

sts_initial.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=nD, r=r, T=T_i_initial)

sts_initial.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)
sts_initial.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_DISABLED)

sts_initial.eqsys.n_i.setJET_CWrecycling()

sts_initial.radialgrid.setB0(B)
sts_initial.radialgrid.setMinorRadius(a)
sts_initial.radialgrid.setMajorRadius(r_0)
sts_initial.radialgrid.setWallRadius(r_wall)
sts_initial.radialgrid.setVesselVolume(V_vessel)
#sts_initial.radialgrid.setElongation(kappa)
sts_initial.radialgrid.setRecyclingCoefficient1(1)
sts_initial.radialgrid.setRecyclingCoefficient2(0)
sts_initial.radialgrid.setRecyclingCoefficient3(1)

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

sts_initial.save('STREAMSettings_initial.h5')



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

sto_list.append(sto_final)

sts_initial.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
sts_initial.eqsys.E_field.setInitialProfile(efield=E_initial)
sts_initial.eqsys.E_field.setInductances(Lp=Lp, Lwall=9.1e-6, M=2.49e-6, Rwall=1e6)
sts_initial.eqsys.E_field.setCircuitVloop(Vloop)

sts_initial.eqsys.T_cold.setType(ColdElectronTemperature.TYPE_SELFCONSISTENT)
sts_initial.eqsys.T_cold.setInitialProfile(T_e_initial)
sts_initial.eqsys.T_cold.setRecombinationRadiation(recombination=Tcold.RECOMBINATION_RADIATION_INCLUDED)

sts_initial.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=nD*(1-a_T), r=r, T=T_i_initial)
sts_initial.eqsys.n_i.addIon(name='T', Z=1, iontype=Ions.IONS_DYNAMIC, n=nD*a_T, r=np.array([0]), tritium=True, T=T_i_initial)

sts_initial.eqsys.n_re.setAvalanche(avalanche=Runaways.AVALANCHE_MODE_NEGLECT)
sts_initial.eqsys.n_re.setDreicer(Runaways.DREICER_RATE_DISABLED)

sts_initial.eqsys.n_i.setJET_CWrecycling()

sts_initial.radialgrid.setB0(B)
sts_initial.radialgrid.setMinorRadius(a)
sts_initial.radialgrid.setMajorRadius(r_0)
sts_initial.radialgrid.setWallRadius(r_wall)
sts_initial.radialgrid.setVesselVolume(V_vessel)
#sts_initial.radialgrid.setElongation(kappa)
sts_initial.radialgrid.setRecyclingCoefficient1(1)
sts_initial.radialgrid.setRecyclingCoefficient2(0)
sts_initial.radialgrid.setRecyclingCoefficient3(1)

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

sts_initial.save('STREAMSettings_initial.h5')



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

sto_list.append(sto_final)

for sto in sto_list:
    V_p = sto.other.stream.V_p[:, 0]

    P_net = P_tot + sto.other.fluid.Tcold_ohmic[:, 0] * V_p

    t = sto.grid.t[:]
    plt.plot(t[1:], P_net, label='Net electron heating power', color='k',linestyle='-')
plt.xlabel('Time [s]')
plt.ylabel('P_net')
plt.legend('D','DT')
plt.show()

for sto in sto_list:
    V_p = sto.other.stream.V_p[:, 0]

    Prad = sto.other.fluid.Tcold_radiation[:, 0] * V_p
    Pequi = sto.other.fluid.Tcold_ion_coll[:, 0] * V_p
    Ptransp = sto.other.scalar.energyloss_T_cold[:, 0] * V_p
    P_tot = Prad + Pequi + Ptransp
    P_net = P_tot + sto.other.fluid.Tcold_ohmic[:, 0] * V_p

    t = sto.grid.t[:]
    plt.plot(t[1:], P_tot, color='r', linestyle='--')
plt.xlabel('Time [s]')
plt.ylabel('P_tot')
plt.legend('D','DT')
plt.show()

for sto in sto_list:
    V_p = sto.other.stream.V_p[:, 0]

    Prad = sto.other.fluid.Tcold_radiation[:, 0] * V_p

    t = sto.grid.t[:]
    plt.plot(t[1:], Prad,  label='Radiation + ionization', color='b', linestyle='-.')
plt.xlabel('Time [s]')
plt.ylabel('P_rad')
plt.legend('D','DT')
plt.show()

for sto in sto_list:
    V_p = sto.other.stream.V_p[:, 0]

    Pequi = sto.other.fluid.Tcold_ion_coll[:, 0] * V_p

    t = sto.grid.t[:]
    plt.plot(t[1:], Pequi, label='Equilibration', color='g', linestyle='--')
plt.xlabel('Time [s]')
plt.ylabel('P_equi')
plt.legend('D','DT')
plt.show()

for sto in sto_list:
    V_p = sto.other.stream.V_p[:, 0]

    Ptransp = sto.other.scalar.energyloss_T_cold[:, 0] * V_p

    t = sto.grid.t[:]
    plt.plot(t[1:], Ptransp, label='Electron transport', color='m', linestyle=':')
plt.xlabel('Time [s]')
plt.ylabel('P_transp')
plt.legend('D','DT')
plt.show()

'''
t = sto_final.grid.t[:]
legend = ['D', 'T']
for sto in sto_list:
    plt.plot(t[1:], P_tot, color='r', linestyle='--')

plt.plot(t[1:], Prad, label='Radiation + ionization', color='b', linestyle='-.')
plt.plot(t[1:], Pequi, label='Equilibration', color='g', linestyle='--')
plt.plot(t[1:], Ptransp, label='Electron transport', color='m', linestyle=':')
plt.plot(t[1:], P_net, label='Net electron heating power', color='k', linestyle='-')
plt.xlabel('Time [s]')
plt.ylabel('Power balance')
#plt.ylim(0, 2e5)
plt.legend()
plt.show()


rad = sto_final.other.fluid.Tcold_radiation[:]
ohmic = sto_final.other.fluid.Tcold_ohmic[:]
transport = sto_final.other.scalar.energyloss_T_cold[:]
equilibration = sto_final.other.fluid.Tcold_ion_coll[:]
plt.plot(sto_final.grid.t[1:], rad+ohmic+transport+equilibration)

#plt.ylim(0,12e5)
plt.xlabel('Time [s]')
plt.ylabel('Total radiation power loss [A]')
plt.show()
'''
# Plasma current
for sto in sto_list:
    sto.eqsys.I_p.plot()

#plt.ylim(0,12e5)
plt.xlabel('Time [s]')
plt.ylabel('Plasma current [A]')
plt.legend('D','DT')
plt.show()



'''
# Total radiation power loss electron
plt.plot(sto_final.grid.t[1:], sto_final.other.stream.V_p[:]*np.transpose(np.array(np.diff(sto_final.eqsys.W_cold[:,0]) / np.diff(sto_final.grid.t[:]))[np.newaxis]))

#plt.ylim(0,12e5)
plt.xlabel('Time [s]')
plt.ylabel('Total radiation power loss [A]')
plt.show()

rad = sto_final.other.fluid.Tcold_radiation
ohmic = sto_final.other.fluid.Tcold_ohmic
transport = sto_final.other.scalar.energyloss_T_cold
equilibration = sto_final.other.fluid.Tcold_ion_coll
plt.plot(sto_final.grid.t[1:], rad+ohmic+transport+equilibration)

#plt.ylim(0,12e5)
plt.xlabel('Time [s]')
plt.ylabel('Total radiation power loss [A]')
plt.show()
'''
'''
# Total radiation power loss ion
plt.plot(sto_final.grid.t[1:], sto_final.other.stream.V_p[:]*np.transpose(np.array(np.diff(sto_final.eqsys.W_i[:,0]) / np.diff(sto_final.grid.t[:]))[np.newaxis]))

#plt.ylim(0,12e5)
plt.xlabel('Time [s]')
plt.ylabel('Total radiation power loss [A]')
plt.show()
'''

# Electron temperature
for sto in sto_list:
    sto_final.eqsys.T_cold.plot()
#plt.ylim(0,12e5)
plt.xlabel('Time [s]')
plt.ylabel('Electron temperature [eV]')
plt.legend('D','DT')
plt.show()


# Electron density
for sto in sto_list:
    sto_final.eqsys.n_cold.plot()
#plt.ylim(0,12e5)
plt.xlabel('Time [s]')
plt.ylabel('Electron density [m^-3]')
plt.legend('D','DT')
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

n_T_0=sto_final.eqsys.n_i['T'][0][:]
legend.append('n_T_0')
plt.plot(sto_final.grid.t[:],n_T_0,'r-')
n_T_1=sto_final.eqsys.n_i['T'][1][:]
legend.append('n_T_1')
plt.plot(sto_final.grid.t[:],n_T_1,'r--')

plt.xlabel('Time [s]')
plt.ylabel('Ion density [m$^{-3}$]')
fontP = FontProperties()
fontP.set_size('x-small')
plt.legend(legend,prop=fontP)
#plt.xlim(0,0.12)
plt.show()

'''