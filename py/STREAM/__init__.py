
try:
    from DREAM import GeriMap
except ModuleNotFoundError:
    # Add DREAM from 
    import sys
    import pathlib
    cwd = pathlib.Path(__file__).resolve()
    p = (cwd.parent / 'extern' / 'DREAM' / 'py').resolve()
    sys.path.append(p)

from . STREAMException import STREAMException
from . STREAMSettings import STREAMSettings
from . STREAMOutput import STREAMOutput
from . runiface import runiface

# Register essential color map
#GeriMap.register()

