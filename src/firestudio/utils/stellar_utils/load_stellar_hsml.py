import numpy as np
import ctypes
import os

def checklen(x): return np.array(x,ndmin=1).shape[0]

def farr(x): return np.array(x,dtype='f',ndmin=1)

def f_cpointer(x): return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

def get_particle_hsml( x, y, z, DesNgb=32, Hmax=0.):
    ## ensure that x,y,z inputs are formatted as 1d float arrays
    x,y,z = farr(x),farr(y),farr(z)

    num_points = checklen(x)

    ## filter out any nans or infs
    ok_mask = np.all(np.isfinite([x,y,z]),axis=0)
    x,y,z = x[ok_mask],y[ok_mask],z[ok_mask]

    if(Hmax==0.):
        max_dx=np.max([
            np.max(x)-np.min(x),
            np.max(y)-np.min(y),
            np.max(z)-np.min(z)]) 

        Hmax=5.*max_dx*(np.float(num_points)**(-1./3.)) ## mean inter-particle spacing

    ## load the routine we need
    ## find the path of the parent directory, removing any symbolic links for safety
    curpath = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    exec_call = os.path.join(curpath,'C_routines/StellarHsml/starhsml.so')
    h_routine = ctypes.cdll[exec_call]

    ## create a C array as output
    H_OUT=(ctypes.c_float*num_points)()

    ## main call to the hsml-finding routine
    h_routine.stellarhsml( 
        ctypes.c_int(num_points), 
        f_cpointer(x),f_cpointer(y),f_cpointer(z),
        ctypes.c_int(DesNgb), 
        ctypes.c_float(Hmax),
        ctypes.byref(H_OUT))

    ## return put the output as a numpy array
    return np.ctypeslib.as_array(np.copy(H_OUT))
