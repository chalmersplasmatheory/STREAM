
try:
    from DREAM.GeriMap import register
except ModuleNotFoundError:
    # Add DREAM from 
    import sys
    import pathlib
    cwd = pathlib.Path('.').resolve()
    p = (cwd.parent / 'extern' / 'DREAM' / 'py').resolve()
    sys.path.append(p)

from . STREAMException import STREAMException
from . STREAMSettings import STREAMSettings
from . STREAMOutput import STREAMOutput

# Register essential color map
GeriMap.register()

