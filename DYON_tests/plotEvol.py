#!/usr/bin/env python3

import matplotlib.pyplot as plt
import os
import sys

sys.path.append('../py')

from STREAM import STREAMOutput


if len(sys.argv) != 2:
    raise Exception('Expected name of parameter to plot as argument.')

param = sys.argv[1]

PATH = 'debug/debugout_1_{}.h5'
iters = []
vals  = []
i = 1
while os.path.exists(PATH.format(i)):
    so = STREAMOutput(PATH.format(i))

    iters.append(i)

    if param == 'n_i':
        vals.append(so.eqsys['n_i']['D'][0][0,0])
    elif param == 'lambda_i':
        vals.append(so.eqsys['lambda_i']['O'][0,0])
    else:
        vals.append(so.eqsys[param][0])

    i += 1


plt.semilogy(iters, vals)
plt.show()

