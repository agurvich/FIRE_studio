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

def snap_ext(snum,four_char=0):
	ext='00'+str(snum);
	if (snum>=10): ext='0'+str(snum)
	if (snum>=100): ext=str(snum)
	if (four_char==1): ext='0'+ext
	if (snum>=1000): ext=str(snum)
	return ext;


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



def load_allstars_hsml(snapdir,snapnum,cosmo=0,use_rundir=0,four_char=0,use_h0=1):
    import gadget

    ## do we already have a file of these pre-calculated?
    rootdir='/n/scratch2/hernquist_lab/phopkins/stellar_hsml_saved/'

    exts=snap_ext(snapnum,four_char=four_char);
    s0=snapdir.split("/"); ss=s0[len(s0)-1]; 
    if(len(ss)==0): ss=s0[len(s0)-2];
    hsmlfile_r=rootdir+ss+'_allstars_hsml_'+exts
    if (use_rundir==1): hsmlfile_r=snapdir+'/allstars_hsml_'+exts
    hsmlfile=hsmlfile_r+'.dat' ## check if binary file exists

    if os.path.exists(hsmlfile): ## it exists! 
        lut=open(hsmlfile,'r'); 
        int_in=array.array('i'); int_in.fromfile(lut,1); nstars=int_in[0];
        h_in=array.array('f'); h_in.fromfile(lut,nstars);         
        lut.close();
        return np.array(np.copy(h_in));
    else: ## no pre-computed file, need to do it ourselves
        have=0; ptype=4;
        ppp=gadget.readsnap(snapdir,snapnum,ptype,h0=use_h0,cosmological=cosmo);
        if(ppp['k']==1): 
            if(ppp['m'].shape[0]>1):
                have=1; pos=ppp['p']; x=pos[:,0]; y=pos[:,1]; z=pos[:,2];
        if (cosmo==0):
            for ptype in [2,3]:
                ppp=gadget.readsnap(snapdir,snapnum,ptype,h0=use_h0,cosmological=0);
                if(ppp['k']==1): 
                    if(ppp['m'].shape[0]>1):
                        pos=ppp['p']
                        if (have==1):
                            x=np.concatenate((x,pos[:,0]));
                            y=np.concatenate((y,pos[:,1]));
                            z=np.concatenate((z,pos[:,2]));
                        else:
                            x=pos[:,0]; y=pos[:,1]; z=pos[:,2];
                        have=1;
        ## ok now we have the compiled positions
        if (have==1):
            h = get_particle_hsml(x,y,z,DesNgb=62);
            ## great now we've got the stars, lets write this to a file for next time
            nstars = checklen(h);
            print nstars
            if (nstars>1):
                lut=open(hsmlfile,'wb');
                lut.write(struct.pack('i',nstars));
                lut.write(struct.pack('f'*len(h),*h));
                lut.close();
                return h;
            else:
                return 0;
        else:
            return 0; 
    return 0; ## failed to find stars

