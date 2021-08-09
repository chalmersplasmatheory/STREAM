import numpy as np
import sys
sys.path.append('../py')
from STREAM import STREAMOutput
import matplotlib.pyplot as plt

do = STREAMOutput('output_final.h5')
Z0=0
plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_posRecombination[:,Z0,0], label='posRec')
plt.plot(do.grid.t[1:],do.other.stream.neutralinflux['D'][:,0]/do.other.stream.V_n_tot['D'][:,0], label='posRec')


#for Z0 in range(9,26):
#    if np.sum(do.other.stream.ionrateequation_posIonization[:,Z0,0])==0:
#        continue
'''
Z0=0
plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_posIonization[:,Z0,0], label='posIoniz')
plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_negIonization[:,Z0,0], label='negIoniz')
plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_posRecombination[:,Z0,0], label='posRec')
plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_negRecombination[:,Z0,0], label='negRec')
plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_posChargeExchange[:,Z0,0], label='posCX')
plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_negChargeExchange[:,Z0,0], label='negCX')
plt.legend()
#'''
plt.show()