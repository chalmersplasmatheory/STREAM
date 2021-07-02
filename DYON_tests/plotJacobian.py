#!/usr/bin/env python3

import os
import sys

sys.path.append(os.environ['DREAMPATH'] + '/tools')
from dreamdebug import *


d = DREAMEqsys('eqsys.txt')

J1 = load('petsc_jac')
J2 = load('petsc_jac_num')


#cmp(J1, J2, eqsys=d, tollow=0.99)
plotrowl(J1, J2, eqsys=d, row=23, legend=['Analytical', 'Numerical'])
#plotrowl(J1, eqsys=d, row=10)

plt.show()
