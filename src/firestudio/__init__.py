import os
import subprocess

def build_C_routines(prepend=''):
    ## cd to whichever folder contains this file
    os.chdir(os.path.dirname(__file__))
    C_routine_subdir = 'utils/C_routines'

    ## build each of the C sub-routines
    C_routines = os.listdir(C_routine_subdir)
    for C_routine in C_routines:
        process = subprocess.Popen(
            "make",
            shell=True,
            cwd=os.path.join(C_routine_subdir,C_routine))
        process.wait()
