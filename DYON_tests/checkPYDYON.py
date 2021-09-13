#!/usr/bin/env python3

import numpy as np
import sys

sys.path.append('../py')

from PYDYON import compareToSTREAM, compareToSTREAMdt
from STREAM import STREAMSettings, STREAMOutput

OUTFILE = 'final'

ss = STREAMSettings(f'settings_{OUTFILE}.h5')
so = STREAMOutput(f'output_{OUTFILE}.h5')

compareToSTREAM(ss, so)
#compareToSTREAMdt(ss, so)

