#!/usr/bin/env python3

import numpy as np
import sys

sys.path.append('../../../py')

from PYDYON import compareToSTREAM, compareToSTREAMdt
from STREAM import STREAMSettings, STREAMOutput

OUTFILE = '2'

ss = STREAMSettings(f'settings{OUTFILE}.h5')
so = STREAMOutput(f'output{OUTFILE}.h5')

compareToSTREAM(ss, so)
#compareToSTREAMdt(ss, so)

