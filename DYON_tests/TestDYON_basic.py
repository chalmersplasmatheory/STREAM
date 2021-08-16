# Bestäm

#'''
import numpy as np
import sys
import matplotlib.pyplot as plt
from PlasmaParameters import evaluateSpitzerConductivity

sys.path.append('../extern/DREAM/py/')
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ColdElectrons as ColdElectrons
import DREAM.Settings.Equations.ColdElectronTemperature as ColdElectronTemperature
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.ElectricField as ElectricField_d
import DREAM.Settings.Equations.ElectricField as ElectricField

sys.path.append('../py/')
from STREAM.STREAMSettings import STREAMSettings
from STREAM.STREAMOutput import STREAMOutput
from STREAM import runiface
import STREAM.Settings.Equations.IonSpecies as Ions
#import STREAM.Settings.Equations.ElectricField as ElectricField
import STREAM.Settings.TransportSettings as Transport
#'''

# Grid parameters

tMax_initial = 1e-4  # simulation time in seconds
Nt_initial   = 4000   # number of time steps
tMax_final   = 1e-1  # simulation time in seconds
Nt_final     = 1000   # number of time steps

pgp_list = np.linspace(5,7,2)#[5e-5, 7e-5]
sto_list = []
#i_list = [5, 7]

for pgp in pgp_list:
    gamma = 2e-3
    n_D_0 = 3.22e22 * pgp *1e-5
    n_D = n_D_0 * np.array([[1-gamma],[gamma]])

    V_vessel = 100
    B        = 2.3
    a        = 0.5 # Instead of EFIT-data 0.08519 (think there's a typo should be 0.8519)
    r_0      = 3   # Instead of EFIT-data 3.0381 (Should use 2.96?)
    r_wall   = 1  # approx 1.5-2.0 m

    c1       = 1
    c2       = 0
    c3       = 1

    r=np.array([0])


    V_loop_wall = 20

    wall_time = 1 / 1e10

    T_e_initial = 1 # eV
    T_i_initial = 0.03
    sigma = evaluateSpitzerConductivity(n=n_D[1], T=T_e_initial, Z=1)
    J = 405.8
    E_initial = J/sigma #V_loop_wall/(2*np.pi*r_0) # Variera från 0 till V_d[0]/(2*np.pi*r_0) och se om simulering är känsligt för denna

    sts_initial = STREAMSettings()

    sts_initial.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
    sts_initial.eqsys.E_field.setInitialProfile(efield=E_initial)
    sts_initial.eqsys.E_field.setBoundaryCondition(ElectricField.BC_TYPE_TRANSFORMER, V_loop_wall_R0=V_loop_wall/r_0, times=0, inverse_wall_time=1/wall_time, R0=r_0)
    #sts_initial.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
    #sts_initial.eqsys.E_field.setInitialProfile(efield=E_initial)
    #sts_initial.eqsys.E_field.setInductances(Lp=8e-6, Lwall=9.1e-6, M=2.49e-6, Rwall=1e6)
    #sts_initial.eqsys.E_field.setCircuitVloop(V_loop_wall, 0)

    sts_initial.eqsys.T_cold.setType(ColdElectronTemperature.TYPE_SELFCONSISTENT)
    sts_initial.eqsys.T_cold.setInitialProfile(T_e_initial)
    sts_initial.eqsys.T_cold.setRecombinationRadiation(recombination=ColdElectronTemperature.RECOMBINATION_RADIATION_INCLUDED)

    sts_initial.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=n_D, r=np.array([0]), T=T_i_initial)

    sts_initial.eqsys.n_re.setAvalanche(False)
    sts_initial.eqsys.n_re.setDreicer(False)

    sts_initial.eqsys.n_i.setJET_CWrecycling()

    sts_initial.radialgrid.setB0(B)
    sts_initial.radialgrid.setMinorRadius(a)
    sts_initial.radialgrid.setMajorRadius(r_0)
    sts_initial.radialgrid.setWallRadius(r_wall)
    sts_initial.radialgrid.setVesselVolume(V_vessel)
    sts_initial.radialgrid.setRecyclingCoefficient1(c1)
    sts_initial.radialgrid.setRecyclingCoefficient2(c2)
    sts_initial.radialgrid.setRecyclingCoefficient3(c3)

    sts_initial.hottailgrid.setEnabled(False)
    sts_initial.runawaygrid.setEnabled(False)

    sts_initial.solver.setType(Solver.NONLINEAR)

    sts_initial.timestep.setTmax(tMax_initial)
    sts_initial.timestep.setNt(Nt_initial)


    sts_initial.solver.preconditioner.setEnabled(False)

    #sts_initial.solver.setVerbose(True)

    sts_initial.other.include('fluid', 'stream')

    #sts_initial.save('STREAMSettings_initial.h5')

    #sts_initial.radialgrid.setWallRadius(r_wall)

    sto_initial = runiface(sts_initial, 'output_initial_STREAM_'+str(pgp)+'.h5',
                           quiet=False)
    #'''
    sts_final = STREAMSettings(sts_initial)
    sts_final.timestep.setTmax(tMax_final)
    sts_final.timestep.setNt(Nt_final)

    sts_final.fromOutput('output_initial_STREAM_'+str(pgp)+'.h5')
    sts_final.output.setFilename('output_final_STREAM_'+str(pgp)+'.h5')

    sto_final = runiface(sts_final, 'output_final_STREAM_'+str(pgp)+'.h5',
                            quiet=False)
    #'''
    sto_list.append(sto_final)
'''
for pgp in pgp_list:
    gamma = 2e-3
    n_D_0 = 3.22e22 * pgp *1e-5
    n_D = n_D_0 * np.array([[1-gamma],[gamma]])

    V_vessel = 100
    B        = 2.3
    a        = 0.5 # Instead of EFIT-data 0.08519 (think there's a typo should be 0.8519)
    r_0      = 3   # Instead of EFIT-data 3.0381 (Should use 2.96?)
    r_wall   = 1  # approx 1.5-2.0 m

    c1       = 1
    c2       = 0
    c3       = 1

    r=np.array([0])


    V_loop_wall = 20

    T_e_initial = 1 # eV
    T_i_initial = 0.03
    sigma = evaluateSpitzerConductivity(n=n_D[1], T=T_e_initial, Z=1)
    J = 405.8
    E_initial = J/sigma #V_loop_wall/(2*np.pi*r_0) # Variera från 0 till V_d[0]/(2*np.pi*r_0) och se om simulering är känsligt för denna

    sts_initial = STREAMSettings()

    #sts_initial.eqsys.E_field.setType(ElectricField.TYPE_SELFCONSISTENT)
    #sts_initial.eqsys.E_field.setInitialProfile(efield=E_initial)
    #sts_initial.eqsys.E_field.setBoundaryCondition(ElectricField.BC_TYPE_TRANSFORMER, V_loop_wall_R0=V_loop_wall/r_0, times=0, inverse_wall_time=1/wall_time, R0=r_0)
    sts_initial.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
    sts_initial.eqsys.E_field.setInitialProfile(efield=E_initial)
    sts_initial.eqsys.E_field.setInductances(Lp=8e-6, Lwall=9.1e-6, M=2.49e-6, Rwall=1e6)
    sts_initial.eqsys.E_field.setCircuitVloop(V_loop_wall, 0)

    sts_initial.eqsys.T_cold.setType(ColdElectronTemperature.TYPE_SELFCONSISTENT)
    sts_initial.eqsys.T_cold.setInitialProfile(T_e_initial)
    sts_initial.eqsys.T_cold.setRecombinationRadiation(recombination=ColdElectronTemperature.RECOMBINATION_RADIATION_INCLUDED)

    sts_initial.eqsys.n_i.addIon(name='D', Z=1, iontype=Ions.IONS_DYNAMIC, n=n_D, r=np.array([0]), T=T_i_initial)

    sts_initial.eqsys.n_re.setAvalanche(False)
    sts_initial.eqsys.n_re.setDreicer(False)

    sts_initial.eqsys.n_i.setJET_CWrecycling()

    sts_initial.radialgrid.setB0(B)
    sts_initial.radialgrid.setMinorRadius(a)
    sts_initial.radialgrid.setMajorRadius(r_0)
    sts_initial.radialgrid.setWallRadius(r_wall)
    sts_initial.radialgrid.setVesselVolume(V_vessel)
    sts_initial.radialgrid.setRecyclingCoefficient1(c1)
    sts_initial.radialgrid.setRecyclingCoefficient2(c2)
    sts_initial.radialgrid.setRecyclingCoefficient3(c3)

    sts_initial.hottailgrid.setEnabled(False)
    sts_initial.runawaygrid.setEnabled(False)

    sts_initial.solver.setType(Solver.NONLINEAR)

    sts_initial.timestep.setTmax(tMax_initial)
    sts_initial.timestep.setNt(Nt_initial)


    sts_initial.solver.preconditioner.setEnabled(False)

    #sts_initial.solver.setVerbose(True)

    sts_initial.other.include('fluid', 'stream')

    #sts_initial.save('STREAMSettings_initial.h5')

    #sts_initial.radialgrid.setWallRadius(r_wall)

    sto_initial = runiface(sts_initial, 'STOs/output_i'+str(pgp)+'.h5',
                           quiet=False)

    sts_final = STREAMSettings(sts_initial)
    sts_final.timestep.setTmax(tMax_final)
    sts_final.timestep.setNt(Nt_final)

    sts_final.fromOutput('STOs/output_i'+str(pgp)+'.h5')
    sts_final.output.setFilename('STOs/output_f'+str(pgp)+'.h5')
    
    sto_final = runiface(sts_final, 'STOs/output_f'+str(pgp)+'.h5',
                            quiet=False)

    sto_list.append(sto_final)
'''
# Power consumption
legend = []
for sto, i in zip(sto_list,pgp_list):
    t = sto.grid.t[:] 
    plt.plot(t[1:], np.diff(sto.eqsys.W_cold[:,0]) / np.diff(sto.grid.t[:]))
    legend.append('PGP='+str(i))
plt.xlabel('Time [s]')
plt.ylabel('Power Consumption [W/m^3]')
plt.legend([])
plt.legend(legend)
#plt.legend(['STREAM 5','STREAM 7', 'DYON 5', 'DYON 7'])
plt.show()


# Plasma current
for sto in sto_list:
    sto.eqsys.I_p.plot()
plt.xlabel('Time [s]')
plt.ylabel('Plasma current [A]')
plt.yscale('log')
plt.legend(legend)
#plt.legend(['STREAM 5','STREAM 7', 'DYON 5', 'DYON 7'])
plt.show()
#'''
