import numpy as np
import ctypes
import math
import os.path
import struct
import array

def checklen(x):
    return len(np.array(x,ndmin=1));
def fcor(x):
    return np.array(x,dtype='f',ndmin=1)
def vfloat(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float));

def ok_scan(input,xmax=1.0e10,pos=0):
    if (pos==1):
        return (np.isnan(input)==False) & (abs(input)<=xmax) & (input > 0.);
    else:
        return (np.isnan(input)==False) & (abs(input)<=xmax);

def get_particle_hsml( x, y, z, DesNgb=32, Hmax=0.):
    x=fcor(x); y=fcor(y); z=fcor(z); N=checklen(x); 
    ok=(ok_scan(x) & ok_scan(y) & ok_scan(z)); x=x[ok]; y=y[ok]; z=z[ok];
    if(Hmax==0.):
        dx=np.max(x)-np.min(x); dy=np.max(y)-np.min(y); dz=np.max(z)-np.min(z); ddx=np.max([dx,dy,dz]); 
        Hmax=5.*ddx*(np.float(N)**(-1./3.)); ## mean inter-particle spacing

    ## load the routine we need
    curpath = os.path.realpath(__file__)
    curpath = curpath[:len("utils")+curpath.index("utils")] #split off this filename
    exec_call=os.path.join(curpath,'stellar_utils/c_libraries/StellarHsml/starhsml.so')
    h_routine=ctypes.cdll[exec_call];

    h_out_cast=ctypes.c_float*N; H_OUT=h_out_cast();
    ## main call to the hsml-finding routine
    h_routine.stellarhsml( ctypes.c_int(N), \
        vfloat(x), vfloat(y), vfloat(z), ctypes.c_int(DesNgb), \
        ctypes.c_float(Hmax), ctypes.byref(H_OUT) )

    ## now put the output arrays into a useful format 
    h = np.ctypeslib.as_array(np.copy(H_OUT));
    return h;
