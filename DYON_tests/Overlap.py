import numpy as np
import sys
sys.path.append('../py')
from STREAM import STREAMOutput
import matplotlib.pyplot as plt

sto_initial = STREAMOutput('output_initial.h5')
sto_final = STREAMOutput('output_final.h5')

plt.plot(sto_initial.grid.t[:],sto_initial.eqsys.E_field[:], label='initial')
plt.plot(sto_final.grid.t[:]+1e-2,sto_final.eqsys.E_field[:], label='Final')
plt.legend()
plt.xlabel('t [s]')
plt.ylabel('E_field [V/m]')
plt.show()

plt.plot(sto_initial.grid.t[:],sto_initial.eqsys.T_cold[:], label='initial')
plt.plot(sto_final.grid.t[:]+1e-2,sto_final.eqsys.T_cold[:], label='Final')
plt.legend()
plt.xlabel('t [s]')
plt.ylabel('T_cold [eV]')
plt.show()

plt.plot(sto_initial.grid.t[:],sto_initial.eqsys.I_p[:], label='initial')
plt.plot(sto_final.grid.t[:]+1e-2,sto_final.eqsys.I_p[:], label='Final')
plt.legend()
plt.xlabel('t [s]')
plt.ylabel('I_p [A]')
plt.show()

plt.plot(sto_initial.grid.t[:],sto_initial.eqsys.I_wall[:], label='initial')
plt.plot(sto_final.grid.t[:]+1e-2,sto_final.eqsys.I_wall[:], label='Final')
plt.legend()
plt.xlabel('t [s]')
plt.ylabel('I_wall [A]')
plt.show()

plt.plot(sto_initial.grid.t[:],sto_initial.eqsys.n_cold[:], label='initial')
plt.plot(sto_final.grid.t[:]+1e-2,sto_final.eqsys.n_cold[:], label='Final')
plt.legend()
plt.xlabel('t [s]')
plt.ylabel('n_cold [m^-3]')
plt.show()

plt.plot(sto_initial.grid.t[:],sto_initial.eqsys.V_loop_w[:], label='initial')
plt.plot(sto_final.grid.t[:]+1e-2,sto_final.eqsys.V_loop_w[:], label='Final')
plt.legend()
plt.xlabel('t [s]')
plt.ylabel('V_loop_w [V]')
plt.show()

'''
plt.plot(sto_initial.grid.t[:],sto_initial.eqsys.n_i['D'][:], label='initial')
plt.plot(sto_final.grid.t[:]+1e-2,sto_final.eqsys.n_i['D'][:], label='Final')
plt.legend()
plt.xlabel('t [s]')
plt.ylabel('n_i_D [m^-3]')
plt.show()

plt.plot(sto_initial.grid.t[:],sto_initial.eqsys.n_i['O'][:], label='initial')
plt.plot(sto_final.grid.t[:]+1e-2,sto_final.eqsys.n_i['O'][:], label='Final')
plt.legend()
plt.xlabel('t [s]')
plt.ylabel('n_i_O [m^-3]')
plt.show()
#'''