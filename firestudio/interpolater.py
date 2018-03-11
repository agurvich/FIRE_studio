import itertools 

from all_utils import * 
from gal_utils import Galaxy

import matplotlib.pyplot as plt
import numpy as np 
import os

from multiprocessing.pool import ThreadPool
from multiprocessing import Pool

def fauxrenderPatch(sub_res,ax,
    frame_center,frame_width,
    frame_depth=None,savefig=0,noaxis=0,
    theta=0,phi=0,psi=0,**kwargs):

    indices = extractRectangularVolumeIndices(sub_res['p'],
        frame_center,frame_width,frame_width if frame_depth is None else frame_depth)

    pos = sub_res['p'] - frame_center # want to rotate about frame_center
    pos_rot = rotateEuler(theta,phi,psi,pos) +frame_center # add back the offset post rotation...?

    xs,ys = pos_rot[:,:2][indices].T

    twoDHist(ax,xs,ys,bins=100)

    ax.set_ylim(frame_center[1]-frame_width,frame_center[1]+frame_width)
    ax.set_xlim(frame_center[0]-frame_width,frame_center[0]+frame_width)

    if noaxis:
        ax.axis('off')
    return ax

def twoDHist(ax,xs,ys,bins,weights=None,norm='',cbar=0):
    if norm=='':
        from matplotlib.colors import LogNorm
        norm=LogNorm()
    cmap=plt.get_cmap('afmhot')
    h,xedges,yedges=np.histogram2d(xs,ys,weights=weights,bins=bins)
    ax.imshow(h.T,cmap=cmap,origin='lower',
    norm=norm,extent=[min(xedges),max(xedges),min(yedges),max(yedges)])
    if cbar:
        plt.colorbar()
    return h,xedges,yedges
    
def rotateEuler(theta,phi,psi,pos):
    ## if need to rotate at all really -__-
    if theta==0 and phi==0 and psi==0:
        return pos
    # rotate particles by angle derived from frame number
    pi        = 3.14159265e0
    theta_rad = pi*theta/ 1.8e2
    phi_rad   = pi*phi  / 1.8e2
    psi_rad   = pi*psi  / 1.8e2

    # construct rotation matrix
    rot_matrix      = new_matrix(3,3)
    #print 'theta = ',theta_rad
    #print 'phi   = ',phi_rad
    #print 'psi   = ',psi_rad
    rot_matrix[0][0] =  np.cos(phi_rad)*np.cos(psi_rad)
    rot_matrix[0][1] = -np.cos(phi_rad)*np.sin(psi_rad)
    rot_matrix[0][2] =  np.sin(phi_rad)
    rot_matrix[1][0] =  np.cos(theta_rad)*np.sin(psi_rad) + np.sin(theta_rad)*np.sin(phi_rad)*np.cos(psi_rad)
    rot_matrix[1][1] =  np.cos(theta_rad)*np.cos(psi_rad) - np.sin(theta_rad)*np.sin(phi_rad)*np.sin(psi_rad)
    rot_matrix[1][2] = -np.sin(theta_rad)*np.cos(phi_rad)
    rot_matrix[2][0] =  np.sin(theta_rad)*np.sin(psi_rad) - np.cos(theta_rad)*np.sin(phi_rad)*np.cos(psi_rad)
    rot_matrix[2][1] =  np.sin(theta_rad)*np.cos(psi_rad) - np.cos(theta_rad)*np.sin(phi_rad)*np.sin(psi_rad)
    rot_matrix[2][2] =  np.cos(theta_rad)*np.cos(phi_rad)
    #print rot_matrix

    ## this needs to be done before passing to rotateeuler
    # translate particles so centre = origin
    #pos[:,0] -= (Xmin + (L_x/2.))
    #pos[:,1] -= (Ymin + (L_y/2.))
    #pos[:,2] -= (Zmin + (L_z/2.))

    ## global variable inside real render 
    n_box = pos.shape[0]

    # rotate about each axis with a matrix operation
    pos_rot = np.ndarray((n_box, 3), dtype=np.float32)
    for ipart in range(n_box):
        pos_matrix = new_matrix(3,1)
        pos_matrix[0][0] = pos[ipart,0]
        pos_matrix[1][0] = pos[ipart,1]
        pos_matrix[2][0] = pos[ipart,2]
        rotated = np.matmul(rot_matrix,pos_matrix)
        pos_rot[ipart,0] = rotated[0][0]
        pos_rot[ipart,1] = rotated[1][0]
        pos_rot[ipart,2] = rotated[2][0]
    
    ## occurs outside rotateeuler
    # translate particles back
    #pos[:,0] += (Xmin + (L_x/2.))
    #pos[:,1] += (Ymin + (L_y/2.))
    #pos[:,2] += (Zmin + (L_z/2.))
    
    return pos_rot

def new_matrix(m,n):
    # Create zero matrix
    matrix = [[0 for row in range(n)] for col in range(m)]
    return matrix

###### BEGIN MOVIE PERSPECTIVE+TIME INTERPOLATION CODE
def drawTimeChangingPath(
    frame_width, # fixed frame width, like FoV
    frame_centers,zooms, # position and how far zoomed in am i? 
    thetas, phis, psis, # rotation angles
    nstepss, # how many interpolation frames should there be between each keyframe?
    steps_per_snap,
    start_snap
    ):
    """ if you want to fly around just a single snapshot,
        just pass steps_per_snap very large"""

    ## exclusive prefix sum
    nframe_offsets = [0]
    for nsteps in nstepss[:-1]:
        nframe_offsets+=[nframe_offsets[-1]+nsteps]

    finsnap = start_snap + (nframe_offsets[-1]+nstepss[-1])/steps_per_snap + 1

    ## setup interpolation bounds between keyframes
    dkeyframes = itertools.izip(
        nstepss,
        frame_centers[:-1],frame_centers[1:],
        frame_width/zooms[:-1],frame_width/zooms[1:],
        thetas[:-1],thetas[1:],
        phis[:-1],phis[1:],
        psis[:-1],psis[1:],
        itertools.repeat(steps_per_snap),
        itertools.repeat(start_snap),
        nframe_offsets,
        )

    for dkeyframe in dkeyframes:
        interpolateKeyFrames(*dkeyframe)
            
def interpolateKeyFrames(
    nsteps,
    r0,rf,
    frame_width0,frame_widthf,
    thetamin,thetamax,
    phimin,phimax,
    psimin,psimax,
    steps_per_snap,
    start_snap,
    nframe_offset=0,
    mps=1
    ):
    """ in general nframe_offset will be nkey * nsteps but
        want to let key frame interpolation go at different speeds..."""

    ## the starting snapshot for this key frame interpolation
    start_snap += nframe_offset/steps_per_snap
    nsnaps_this_interpolation = int(np.ceil(1.0*nsteps/steps_per_snap))

    ## figure out how many frames of the interpolation each snapshot should be getting
    ## correct for frames this snapshot has contributed to as part of the previous
    ## key frame interpolation
    frames_to_do = np.array([ 
        steps_per_snap-(nframe_offset%steps_per_snap)*(i==0) for i in xrange(nsnaps_this_interpolation)])
    

    ## clamp the right edge to make sure we don't step too far
    frames_to_do[-1]=nsteps-sum(frames_to_do[:-1])
    
    ## i feel like there's a smarter way to correct for this edge case
    if frames_to_do[-1]>steps_per_snap:
        frames_to_do=np.append(frames_to_do[:-1],[steps_per_snap,frames_to_do[-1]-steps_per_snap])
        nsnaps_this_interpolation+=1

    ## make sure we haven't accidentally done too many frames
    assert sum(frames_to_do) == nsteps

    snapnums = start_snap+np.arange(nsnaps_this_interpolation)
    print frames_to_do
    sub_step_boundss = sub_interpolate_values(0,23,frames_to_do)

    
    xmin,ymin,zmin = r0
    xmax,ymax,zmax = rf

    ## split up the interpolation into bins
    sub_interp_valsss = [
        sub_interpolate_values(lval,rval,frames_to_do) for lval,rval in 
            [(xmin,xmax),
            (ymin,ymax),
            (zmin,zmax),
            (frame_width0,frame_widthf),
            (thetamin,thetamax),
            (phimin,phimax),
            (psimin,psimax),
            ]
        ] 

    ## exclusive prefix sum
    nframe_offsets = [nframe_offset]
    for nsteps in frames_to_do[:-1]:
        nframe_offsets+=[nframe_offsets[-1]+nsteps]


    ## group arguments together to be passed to multiprocess frame drawing
    argss = itertools.izip(
        snapnums, 
        frames_to_do,
        nframe_offsets, ## so we know what to save the frame number as
        *sub_interp_valsss ## two dimensions are unwrapped here, one by zip and one by *
    )

    ## check that my arguments are being split correctly...
    """
    for i,args in enumerate(argss):
        print '------'
        print i,args
        print '------'

    print
    """


    if mps: 
        my_thread_pool = ThreadPool(min(len(snapnums),5)) 
        my_thread_pool.map(multiProcessFrameDrawingWrapper,argss)
        my_thread_pool.close()
    else:
        for args in argss:
            multiProcessFrameDrawingWrapper(args)

def sub_interpolate_values(beg,fin,stepss):
    """helper function to split interpolation up into binned sub-steps"""
    interpolation_steps = np.linspace(beg,fin,np.sum(stepss))
    sub_interpolation_bounds = []
    for i,steps in enumerate(stepss):
        place = sum(stepss[:i])
        relevant_steps = interpolation_steps[place:place+steps]
        sub_interpolation_bounds += [relevant_steps]#[(relevant_steps[0],relevant_steps[-1])]
    return sub_interpolation_bounds

def multiProcessFrameDrawingWrapper(args):
    """wrapper function to pass to multiprocess.pool.map"""
    return multiProcessFrameDrawing(*args)

def multiProcessFrameDrawing(
    snapnum,nsteps_this_snap,
    nframe_offset,
    xs,ys,zs,
    frame_widths,
    thetas,phis,psis):
    """ Loads snapshot and draws 'relevant' (see above) interpolated frames from it"""

    ## repack coordinates into vector
    frame_centers = np.array([xs,ys,zs]).T

    ## load the galaxy
    snapdir = "/home/abg6257/projects/isoDisk/makeNewDisk_tests/rescaled_snonly/rescaled_fiducial/output"
    galaxy = Galaxy(
        'isogal_snonly',
        snapdir,
        snapnum,
        cosmological=0,
        datadir='/projects/b1026/agurvich/data',easy_load=0)

    ## let's assume that we're only loading the gas particles here
    galaxy.load_gas()

    ##  for an isolated galaxy, not doing much
    galaxy.extractMainHalo()

    drawFrames(
        galaxy,frame_centers,frame_widths,
        thetas,phis,psis,mps=1,offset=nframe_offset)
    
def drawFrames(galaxy,frame_centers,frame_widths,thetas,phis,psis,offset=0,mps=0):
    nframes = len(frame_centers)
    args = itertools.izip(
            range(nframes),
            itertools.repeat(galaxy.sub_res),
            frame_centers,frame_widths,
            thetas,phis,psis,
            itertools.repeat(offset))

    if mps: 
        my_pool = Pool(25)
        my_pool.map(multiWrapper,args)
        my_pool.close()
            
    else:
        for arg in args:
            drawFrame(galaxy,*arg)

def multiWrapper(args):
    multidrawFrame(*args)

def multidrawFrame(i,sub_res,frame_center,frame_width,theta,phi,psi,offset):
    """ uses fauxrenderPatch for testing purposes, or full 'render galaxy' 
        if you want to spend the time on it..."""
    ax = plt.gca()
    fauxrenderPatch(sub_res,
        ax,frame_center,frame_width,
        theta=theta,phi=phi,psi=psi,noaxis=1)
    plt.savefig(os.path.join(sub_res['datadir'],'frame_%04d.png'%(i+offset)))
    plt.clf()

def drawFrame(galaxy,i,frame_center,frame_width,theta,phi,psi,offset):
    ax = plt.gca()
    fauxrenderPatch(
        galaxy.sub_res,
        ax,frame_center,frame_width,
        theta=theta,phi=phi,psi=psi,noaxis=1)
    plt.savefig(os.path.join(galaxy.datadir,'frame_%04d.png'%(i+offset)))
    plt.clf()
