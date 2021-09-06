import numpy as np
import sys
sys.path.append('../py')
from STREAM import STREAMOutput
import matplotlib.pyplot as plt

do = STREAMOutput('output_final.h5')

"""
Z0=0
plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_posRecombination[:,Z0,0], label='posRec')
plt.plot(do.grid.t[1:],do.other.stream.neutralinflux['D'][:,0]/do.other.stream.V_n_tot['D'][:,0], label='posRec')
#"""
#for Z0 in range(9,26):
#    if np.sum(do.other.stream.ionrateequation_posIonization[:,Z0,0])==0:
#        continue
for Z0 in range(0,18):
    plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_posIonization[:,Z0,0], label='posIoniz')
    plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_negIonization[:,Z0,0], label='negIoniz')
    plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_posRecombination[:,Z0,0], label='posRec')
    plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_negRecombination[:,Z0,0], label='negRec')
    plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_posChargeExchange[:,Z0,0], label='posCX')
    plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_negChargeExchange[:,Z0,0], label='negCX')
    if Z0 > 8:
        print(str(Z0-9))
    else:
        if Z0 > 1:
            print(str(Z0 - 2))
        else:
            print(str(Z0))
    if Z0 == 0:
        plt.plot(do.grid.t[1:],do.other.stream.neutralinflux['D'][:,0]/do.other.stream.V_n_tot['D'][:,0], label='transport')
    if Z0 == 2:
        plt.plot(do.grid.t[1:], do.other.stream.neutralinflux['C'][:, 0]/ do.other.stream.V_n_tot['C'][:, 0],
                 label='transport')
    if Z0 == 9:
        plt.plot(do.grid.t[1:], do.other.stream.neutralinflux['O'][:, 0] / do.other.stream.V_n_tot['O'][:, 0],
                 label='transport')
    plt.legend()
    
    plt.show()
    dndt = do.other.stream.ionrateequation_posIonization[:,Z0,0]\
         + do.other.stream.ionrateequation_negIonization[:,Z0,0]\
         + do.other.stream.ionrateequation_posRecombination[:,Z0,0]\
         + do.other.stream.ionrateequation_negRecombination[:,Z0,0] \
         + do.other.stream.ionrateequation_posChargeExchange[:,Z0,0]\
         + do.other.stream.ionrateequation_negChargeExchange[:,Z0,0]
    if Z0 == 0:
        dndt+=do.other.stream.neutralinflux['D'][:,0]/do.other.stream.V_n_tot['D'][:,0]
    if Z0 == 2:
        dndt+=do.other.stream.neutralinflux['C'][:, 0]/ do.other.stream.V_n_tot['C'][:, 0]
    if Z0 == 9:
        dndt+=do.other.stream.neutralinflux['O'][:, 0] / do.other.stream.V_n_tot['O'][:, 0]
    plt.plot(do.grid.t[1:],dndt/1e4)
    plt.show()

'''
do = STREAMOutput('STOs/output_DYON_f5.0.h5')
for Z0 in range(0,2):
    plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_posIonization[:,Z0,0], label='posIoniz')
    plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_negIonization[:,Z0,0], label='negIoniz')
    plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_posRecombination[:,Z0,0], label='posRec')
    plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_negRecombination[:,Z0,0], label='negRec')
    #plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_posChargeExchange[:,Z0,0], label='posCX')
    #plt.plot(do.grid.t[1:],do.other.stream.ionrateequation_negChargeExchange[:,Z0,0], label='negCX')
    if Z0 == 0:
        plt.plot(do.grid.t[1:],do.other.stream.neutralinflux['D'][:]/do.other.stream.V_n_tot['D'][:], label='transport')
    plt.legend()

    plt.show()
    dndt = do.other.stream.ionrateequation_posIonization[:,Z0,0]\
         + do.other.stream.ionrateequation_negIonization[:,Z0,0]\
         + do.other.stream.ionrateequation_posRecombination[:,Z0,0]\
         + do.other.stream.ionrateequation_negRecombination[:,Z0,0] \
         + do.other.stream.ionrateequation_posChargeExchange[:,Z0,0]\
         + do.other.stream.ionrateequation_negChargeExchange[:,Z0,0]
    if Z0 == 0:
        dndt+=do.other.stream.neutralinflux['D'][:]/do.other.stream.V_n_tot['D'][:]*3.5e20/1.75e20
    plt.plot(do.grid.t[1:],dndt)
    plt.show()
#'''