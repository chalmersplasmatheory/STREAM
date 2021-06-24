import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.append('../extern/DREAM/py')
import DREAM.Settings.Equations.RunawayElectrons as Runaways
import DREAM.Settings.Solver as Solver
import DREAM.Settings.CollisionHandler as Collisions
import DREAM.Settings.Equations.ColdElectrons as ColdElectrons
import DREAM.Settings.Equations.ColdElectronTemperature as ColdElectronTemperature
import DREAM.Settings.Equations.ElectricField as ElectricField
import DREAM.Settings.Equations.DistributionFunction as DistFunc

sys.path.append('../py/')
from STREAM import runiface
from STREAM.STREAMSettings import STREAMSettings
from STREAM.STREAMOutput import STREAMOutput
import STREAM.Settings.Equations.IonSpecies as Ions
import STREAM.Settings.Solver as Solver
import STREAM.Settings.TransportSettings as Transport

