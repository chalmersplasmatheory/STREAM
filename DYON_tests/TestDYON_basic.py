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


sys.path.append('../py/')
from STREAM.STREAMSettings import STREAMSettings
from STREAM.STREAMOutput import STREAMOutput
from STREAM import runiface
import STREAM.Settings.Equations.IonSpecies as Ions
import STREAM.Settings.Equations.ElectricField as ElectricField
#import STREAM.Settings.Equations.ElectricField as ElectricField
import STREAM.Settings.TransportSettings as Transport
#'''

# Grid parameters

tMax_initial = 1e-4  # simulation time in seconds
Nt_initial   = 4000   # number of time steps
tMax_final   = 5e-1  # simulation time in seconds
Nt_final     = 1e4*tMax_final   # number of time steps

pgp_list = np.linspace(5,7.0,2)
sto_list = []

for pgp in pgp_list:
    gamma = 2e-3
    n_D_0 = 3.22e22 * pgp * 1e-5
    n_D = n_D_0 * np.array([[1-gamma],[gamma]])

    V_vessel = 100
    B        = 2.3
    a        = 0.5
    r_0      = 3
    r_wall   = 1

    c1       = 1
    c2       = 0
    c3       = 1

    r=np.array([0])


    V_loop_wall = 20

    wall_time = 1 / 1e12

    T_e_initial = 1
    T_i_initial = 0.03
    sigma = evaluateSpitzerConductivity(n=n_D[1], T=T_e_initial, Z=1)
    J = 405.8
    E_initial = J/sigma

    sts_initial = STREAMSettings()

    sts_initial.eqsys.E_field.setType(ElectricField_d.TYPE_SELFCONSISTENT)
    sts_initial.eqsys.E_field.setInitialProfile(efield=E_initial)
    sts_initial.eqsys.E_field.setBoundaryCondition(ElectricField_d.BC_TYPE_TRANSFORMER, V_loop_wall_R0=V_loop_wall/r_0, times=0, inverse_wall_time=1/wall_time, R0=r_0)

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

    sts_initial.other.include('fluid', 'stream', 'scalar')

    sto_initial = runiface(sts_initial, 'STOs/output_STREAM_i_'+str(pgp)+'.h5',
                           quiet=False)

    sts_final = STREAMSettings(sts_initial)
    sts_final.timestep.setTmax(tMax_final)
    sts_final.timestep.setNt(Nt_final)
    sts_final.fromOutput('STOs/output_STREAM_i_'+str(pgp)+'.h5')
    sts_final.output.setFilename('STOs/output_STREAM_f_'+str(pgp)+'.h5')
    sto_final = runiface(sts_final, 'STOs/output_STREAM_f_'+str(pgp)+'.h5',
                            quiet=False)

    sto_list.append(sto_final)
#'''
for pgp in pgp_list:
    gamma = 2e-3
    n_D_0 = 3.22e22 * pgp * 1e-5
    n_D = n_D_0 * np.array([[1-gamma],[gamma]])

    V_vessel = 100
    B        = 2.3
    a        = 0.5
    r_0      = 3
    r_wall   = 1

    c1       = 1
    c2       = 0
    c3       = 1

    r=np.array([0])


    V_loop_wall = 20

    T_e_initial = 1
    T_i_initial = 0.03
    sigma = evaluateSpitzerConductivity(n=n_D[1], T=T_e_initial, Z=1)
    J = 405.8
    E_initial = J/sigma

    sts_initial = STREAMSettings()

    sts_initial.eqsys.E_field.setType(ElectricField.TYPE_CIRCUIT)
    sts_initial.eqsys.E_field.setInitialProfile(efield=E_initial)
    sts_initial.eqsys.E_field.setInductances(Lp=8e-6, Lwall=9.1e-6, M=2.49e-6, Rwall=1e6) # Lp=8.65e-6 för likhet med STREAM
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

    sts_initial.other.include('fluid', 'stream', 'scalar')

    sto_initial = runiface(sts_initial, 'STOs/output_DYON_i'+str(pgp)+'.h5',
                           quiet=False)

    sts_final = STREAMSettings(sts_initial)
    sts_final.timestep.setTmax(tMax_final)
    sts_final.timestep.setNt(Nt_final)

    sts_final.fromOutput('STOs/output_DYON_i'+str(pgp)+'.h5')
    sts_final.output.setFilename('STOs/output_DYON_f'+str(pgp)+'.h5')

    sto_final = runiface(sts_final, 'STOs/output_DYON_f'+str(pgp)+'.h5',
                            quiet=False)

    sto_list.append(sto_final)
#'''
# Power consumption
#'''
legend = []
colour=['b--','r--','b-','r-']
for sto, c in zip(sto_list,colour):
    t = sto.grid.t[:]
    plt.plot(t[1:], np.diff(sto.eqsys.W_cold[:,0]) / np.diff(sto.grid.t[:]),c)
    #legend.append('PGP='+str(i)+'e-5 Torr')
plt.xlabel('Time [s]')
plt.ylabel('Power Consumption [W/m^3]')
plt.legend([])
#plt.legend(legend)
plt.legend(['STREAM 5','STREAM 7', 'DYON 5', 'DYON 7'])
plt.show()


# Plasma current
for sto, c in zip(sto_list,colour):
    #sto.eqsys.I_p.plot()
    t = sto.grid.t[:]
    plt.plot(t[:], sto.eqsys.I_p[:], c)
plt.xlabel('Time [s]')
plt.ylabel('Plasma current [A]')
plt.yscale('log')
#plt.legend(legend)
plt.legend(['STREAM 5','STREAM 7', 'DYON 5', 'DYON 7'])
plt.show()


# Degree of ionization
for sto, c in zip(sto_list,colour):
    gamma = sto.eqsys.n_cold[:]/(sto.eqsys.n_cold[:]+sto.eqsys.n_i['D'][0][:])
    plt.plot(sto.grid.t[:],gamma,c)
plt.xlabel('Time [s]')
plt.xlim(0,0.03)
plt.ylabel('Degree of ionization [%]')
#plt.legend(legend)
plt.legend(['STREAM 5','STREAM 7', 'DYON 5', 'DYON 7'])
plt.show()
#'''
#'''
# Electron temperature
for sto, c in zip(sto_list,colour):
    t = sto.grid.t[:]
    #sto.eqsys.T_cold.plot()
    plt.plot(t[:], sto.eqsys.T_cold[:], c)
plt.xlabel('Time [s]')
plt.xlim(0,0.03)
plt.ylabel('Electron temperature [eV]')
#plt.legend(legend)
plt.legend(['STREAM 5','STREAM 7', 'DYON 5', 'DYON 7'])
plt.show()

# Electric field
for sto, c in zip(sto_list,colour):
    t = sto.grid.t[:]
    #sto.eqsys.E_field.plot()
    plt.plot(t[:], sto.eqsys.E_field[:], c)
plt.xlabel('Time [s]')
plt.ylabel('Electric field [V/m]')
#plt.legend(legend)
plt.legend(['STREAM 5','STREAM 7', 'DYON 5', 'DYON 7'])
plt.show()
#'''
# Electron density
for sto, c in zip(sto_list,colour):
    t = sto.grid.t[:]
    #sto.eqsys.n_cold.plot()
    plt.plot(t[:], sto.eqsys.n_cold[:], c)
plt.xlabel('Time [s]')
plt.xlim(0,0.03)
plt.ylabel('Electron density [m^-3]')
#plt.legend(legend)
plt.legend(['STREAM 5','STREAM 7', 'DYON 5', 'DYON 7'])
plt.show()

# Ohmic heating
for sto, c in zip(sto_list, colour):
    t = sto.grid.t[:]
    #sto.other.fluid.Tcold_ohmic.plot()
    V_p = np.transpose(np.array(sto.other.stream.V_p[:, 0])[np.newaxis])
    plt.plot(t[1:], -sto.other.fluid.Tcold_ohmic[:]*V_p,c)
plt.xlabel('Time [s]')
#plt.xlim(0, 0.03)
#plt.ylabel('Ohmic heating [eV]')
# plt.legend(legend)
plt.legend(['STREAM 5','STREAM 7', 'DYON 5', 'DYON 7'])
plt.show()

for sto, c in zip(sto_list, colour):
    t = sto.grid.t[:]
    #sto.other.fluid.Tcold_ohmic.plot()
    V_p = np.transpose(np.array(sto.other.stream.V_p[:, 0])[np.newaxis])
    plt.plot(t[1:], -sto.other.fluid.Tcold_ohmic[:]*V_p,c)
    I_p = sto.eqsys.I_p[1:]
    R_p = 5e-5 * sto.other.fluid.lnLambdaT[:] * sto.other.fluid.Zeff[:] * 2*3/0.5**2 * sto.eqsys.T_cold[1:]**(-3/2)
    P_ohm_DYON = I_p**2*R_p#/V_p
    plt.plot(t[1:],P_ohm_DYON,'k')
    plt.xlabel('Time [s]')
    #plt.xlim(0, 0.05)
    #plt.ylim(0,1.8e5)
    plt.ylabel('Ohmic heating [eV]')
    # plt.legend(legend)
    plt.legend(['STREAM','Avhandling'])
    plt.show()
#plt.xlabel('Time [s]')
#plt.xlim(0,0.03)
#plt.ylabel('Electron temperature [eV]')
#plt.legend(legend)

#plt.show()
'''
# Confinement time
for sto in sto_list:
    sto.other.stream.tau_D.plot()
plt.legend(['STREAM 5','STREAM 7', 'DYON 5', 'DYON 7'])
plt.xlabel('Time [s]')
plt.ylabel('Confinement time [s]')
plt.show()

for sto in sto_list:
    sto.other.stream.tau_D_par.plot()
plt.legend(['STREAM 5','STREAM 7', 'DYON 5', 'DYON 7'])
plt.xlabel('Time [s]')
plt.ylabel('Confinement time [s]')
plt.show()

for sto in sto_list:
    sto.other.stream.tau_D_perp.plot()
plt.legend(['STREAM 5','STREAM 7', 'DYON 5', 'DYON 7'])
plt.xlabel('Time [s]')
plt.ylabel('Confinement time [s]')
plt.show()

# Energy balance
for sto, i in zip(sto_list,pgp_list):
    sto.eqsys.T_cold.plotEnergyBalance(r=0)
    plt.xlabel('Time [s]')
    plt.ylabel('Energy balance for pgp='+str(i)+'e-5 Torr')
    #plt.legend(legend)
    #plt.legend(['STREAM 5','STREAM 7', 'DYON 5', 'DYON 7'])
    plt.show()
#'''

# Power balance
for sto, i in zip(sto_list,pgp_list):
    V_p = sto.other.stream.V_p[:, 0]

    Prad = sto.other.fluid.Tcold_radiation[:, 0] * V_p
    Pequi = sto.other.fluid.Tcold_ion_coll[:, 0] * V_p
    Ptransp = sto.other.scalar.energyloss_T_cold[:, 0] * V_p
    P_tot = Prad + Pequi + Ptransp
    P_net = P_tot + sto.other.fluid.Tcold_ohmic[:, 0] * V_p

    t = sto.grid.t[:]
    plt.plot(t[1:], P_tot, label='Total electron power loss', color='r', linestyle='--')
    plt.plot(t[1:], Prad,  label='Radiation + ionization', color='b', linestyle='-.')
    plt.plot(t[1:], Pequi, label='Equilibration', color='g', linestyle='--')
    plt.plot(t[1:], Ptransp, label='Electron transport', color='m', linestyle=':')
    plt.plot(t[1:], P_net, label='Net electron heating power', color='k',linestyle='-')
    plt.xlabel('Time [s]')
    plt.ylabel('Power balance')
    plt.xlim(0, 0.05)
    plt.ylim(0,2e5)
    #plt.legend()
    plt.show()

# Total amount of deuterium ions
legend = []
colour=['b-','r-','b-','r-']
for sto, c in zip(sto_list, colour):
    n_D_0=np.array(sto.eqsys.n_i['D'][0][1:])
    #legend.append('PGP='+str(i)+', N_D_0')
    plt.plot(sto.grid.t[1:],n_D_0,c)
colour=['b--','r--','b--','r--']
for sto, c in zip(sto_list, colour):
    n_D_1=sto.eqsys.n_i['D'][1][1:]
    plt.plot(sto.grid.t[1:],n_D_1,c)
    #legend.append('PGP='+str(i)+', N_D_1')

plt.xlabel('Time [s]')
plt.ylabel('Total amount of deuterium ions')
plt.legend(['STREAM 5 N_D_0', 'STREAM 7 N_D_0', 'DYON 5 N_D_0', 'DYON 7 N_D_0', 'STREAM 5 N_D_1', 'STREAM 7 N_D_1', 'DYON 5 N_D_1', 'DYON 7 N_D_1'])
plt.legend(legend)
plt.show()
'''
# Total amount of deuterium ions
legend = []
colour=['b-','r-','b--','r--']
for sto, c in zip(sto_list, colour):
    n_D_0=np.array(sto.eqsys.n_i['D'][0][1:])
    n_D_0 = np.array(n_D_0).flatten()
    dnD0dt = np.gradient(n_D_0)
    #legend.append('PGP='+str(i)+', N_D_0')
    plt.plot(sto.grid.t[1:],dnD0dt,c)
plt.plot(np.linspace(0,tMax_final,2),np.linspace(0,0,2),'k')
plt.xlabel('Time [s]')
plt.ylabel('Derivative of ion density [m$^{-3}$]')
plt.legend(['STREAM 5 N_D_0', 'STREAM 7 N_D_0', 'DYON 5 N_D_0', 'DYON 7 N_D_0'])
plt.show()
'''

# Total amount of deuterium ions
legend = []
colour=['b-','r-','b-','r-']
for sto, c in zip(sto_list, colour):
    N_D_0=np.array(sto.eqsys.n_i['D'][0][1:])*np.transpose(np.array(sto.other.stream.V_n_tot['D'][:])[np.newaxis])
    #legend.append('PGP='+str(i)+', N_D_0')
    plt.plot(sto.grid.t[1:],N_D_0,c)
colour=['b--','r--','b--','r--']
for sto, c in zip(sto_list, colour):
    N_D_1=sto.eqsys.n_i['D'][1][1:]*sto.other.stream.V_p[:]
    plt.plot(sto.grid.t[1:],N_D_1,c)
    #legend.append('PGP='+str(i)+', N_D_1')

plt.xlabel('Time [s]')
plt.ylabel('Total amount of deuterium ions')
plt.legend(legend)
plt.legend(['STREAM 5 N_D_0', 'STREAM 7 N_D_0', 'DYON 5 N_D_0', 'DYON 7 N_D_0', 'STREAM 5 N_D_1', 'STREAM 7 N_D_1', 'DYON 5 N_D_1', 'DYON 7 N_D_1'])
plt.show()

# Total amount
legend = []
colour=['b--','r--','b-','r-']
for sto, c in zip(sto_list, colour):
    N_D_0=sto.eqsys.n_i['D'][0][1:]*np.transpose(np.array(sto.other.stream.V_n_tot['D'][:][np.newaxis]))
    N_D_1=sto.eqsys.n_i['D'][1][1:]*sto.other.stream.V_p[:]
    N_D=N_D_0+N_D_1
    plt.plot(sto.grid.t[1:],N_D,c)

plt.xlabel('Time [s]')
plt.ylabel('Total amount of particles')
#fontP = FontProperties()
#fontP.set_size('x-small')
plt.legend(['STREAM 5','STREAM 7', 'DYON 5', 'DYON 7'])
plt.legend(legend)#,prop=fontP)
#plt.xlim(0,0.12)
plt.show()

