# Simple wrapper for running 'streami'

import os
import pathlib
import subprocess
import tempfile

from . STREAMException import STREAMException
from . STREAMOutput import STREAMOutput
from . STREAMSettings import STREAMSettings


STREAMPATH = None

def locatestream():
    global STREAMPATH

    try:
        STREAMPATH = os.environ['STREAMPATH']

        if STREAMPATH[-1] == '/':
            STREAMPATH = STREAMPATH[:-1]

    except KeyError: pass

    if STREAMPATH is None:
        STREAMPATH = (pathlib.Path(__file__).parent / '..' / '..').resolve().absolute()

        if not os.path.isfile('{}/build/iface/streami'.format(STREAMPATH)):
            #raise STREAMException("Unable to locate the STREAMi executable. Try to set the 'STREAMPATH' environment variable.")
            print("WARNING: Unable to locate the STREAMi executable. Try to set the 'STREAMPATH' environment variable.")


def runiface(settings, outfile=None, quiet=False):
    """
    Run 'streami' with the specified settings (which may be either
    a 'STREAMSettings' object or the name of a file containing the
    settings).

    settings: 'STREAMSettings' object or name of file containing settings.
    outfile:  Name of file to write output to (default: 'output.h5')
    """
    global STREAMPATH

    deleteOutput = False
    if outfile is None:
        deleteOutput = True
        outfile = next(tempfile._get_candidate_names())+'.h5'

    infile = None
    if isinstance(settings, STREAMSettings):
        infile = next(tempfile._get_candidate_names())+'.h5'
        settings.output.setFilename(outfile)
        settings.save(infile)
    else:
        infile = settings

    errorOnExit = 0
    p = None
    obj = None
    stderr_data = None
    try:
        if quiet:
            p = subprocess.Popen(['{}/build/iface/streami'.format(STREAMPATH), infile], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        else:
            p = subprocess.Popen(['{}/build/iface/streami'.format(STREAMPATH), infile], stderr=subprocess.PIPE)

        stderr_data = p.communicate()[1].decode('utf-8')

        if p.returncode != 0:
            errorOnExit = 1
        else:
            obj = STREAMOutput(outfile)

            if deleteOutput:
                os.remove(outfile)

    except KeyboardInterrupt:
        errorOnExit = 2
    finally:
        os.remove(infile)

    if errorOnExit == 1:
        print(stderr_data)
        raise STREAMException("STREAMi exited with a non-zero exit code: {}".format(p.returncode))
    elif errorOnExit == 2:
        raise STREAMException("STREAMi simulation was cancelled by the user.")
    else:
        return obj

locatestream()

