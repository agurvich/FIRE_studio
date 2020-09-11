import os
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from abg_python.snapshot_utils import convertSnapToDF

from abg_python.galaxy.gal_utils import Galaxy

from firestudio.studios.gas_studio import GasStudio



def control_flow(times,pair_snapnums,render_kwargs,**galaxy_kwargs):
    prev_galaxy,next_galaxy = None,None
    prev_snapnum,next_snapnum=None,None
    for i,pair in enumerate(pair_snapnums):
        ## keep the current snapnum
        if pair[0] == prev_snapnum:
            pass
        ## step forward in time, swap pointers
        elif pair[0] == next_snapnum:
            prev_galaxy = next_galaxy
            next_galaxy = None
        ## will need to read from disk
        else:
            prev_galaxy = None
        
        if pair[1] == next_snapnum:
            pass
        else:
            next_galaxy = None
            
        ## unpack the pair
        prev_snapnum,next_snapnum = pair
        
        ## load galaxies that are set to None
        prev_galaxy,next_galaxy,changed = load_gals_from_disk(
            prev_snapnum,next_snapnum,
            prev_galaxy,next_galaxy,**galaxy_kwargs)
        
        if changed:
            ## make an interpolated snapshot with these galaxies
            t0,t1 =  prev_galaxy.current_time_Gyr,next_galaxy.current_time_Gyr
            t = times[i]
            time_merged_df = index_match_snapshots_with_dataframes(
                prev_galaxy.sub_snap,
                next_galaxy.sub_snap,
                extra_keys_to_extract=['Temperature'])
            interp_snap = make_interpolated_snap(t,time_merged_df,t0,t1)
            
        print('----- interpolated snap available -----')
        ## let's put the FIREstudio projections into a sub-directory of our Galaxy class instance
        studio_datadir = os.path.join(os.path.dirname(next_galaxy.datadir),'firestudio')
        my_gasStudio = GasStudio(
            studio_datadir,
            next_snapnum,
            next_galaxy.name,
            gas_snapdict=interp_snap,
        )
        
        my_gasStudio.this_setup_id += "_time%.3f"%t
        plt.figure()
        my_gasStudio.render(plt.gca(),**render_kwargs)

def find_bordering_snapnums(snap_times_gyr,dGyr=.005,tmin=None,tmax=None):

    ## handle default maximum time
    tmax = snap_times_gyr[-1] if tmax is None else tmax
    
    ## handle default minimum time
    if tmin is None:
        tmin = snap_times_gyr[0]
    ## remove dGyr so that tmin is included in arange below
    elif tmin - dGyr > 0:
        tmin = tmin-dGyr
    ## create list of times, -1e-9 to avoid landing exactly on a snapshot number
    times_gyr = np.arange(tmax,tmin,-dGyr)[::-1]-1e-9
    
    inds_next = np.argmax((times_gyr - snap_times_gyr[:,None]) < 0 ,axis=0)
    inds_prev = inds_next-1
    return times_gyr,np.array(list(zip(inds_prev,inds_next)))

def load_gals_from_disk(
    prev_snapnum,next_snapnum,
    prev_galaxy,next_galaxy,
    **kwargs):
    
    changed = False
    if prev_galaxy is None:
        #print('loading',prev_snapnum,'from disk')
        prev_galaxy = Galaxy(snapnum=prev_snapnum,**kwargs)
        prev_galaxy.extractMainHalo()
        changed = True
    if next_galaxy is None:
        #print('loading',next_snapnum,'from disk')
        next_galaxy = Galaxy(snapnum=next_snapnum,**kwargs)
        next_galaxy.extractMainHalo()
        changed = True
        
    return prev_galaxy,next_galaxy,changed
        

def linear_interpolate(
    x0,x1,
    t0,t1,
    t):
    return x0 + (x1-x0)/(t1-t0)*(t-t0)

def make_interpolated_snap(t,time_merged_df,t0,t1):
    interp_snap = {}
    for key in time_merged_df.keys():
        if '_next' in key:
            continue
        elif key in ['coord_xs','coord_ys','coord_zs']:
            continue
        interp_snap[key] = linear_interpolate(
            getattr(time_merged_df,key),
            getattr(time_merged_df,key+'_next'),
            t0,t1,
            t).values
    
    ## handle coordinates explicitly
    coords = np.zeros((time_merged_df.shape[0],3))
    for i,key in enumerate(['coord_xs','coord_ys','coord_zs']):
        coords[:,i] = linear_interpolate(
            getattr(time_merged_df,key),
            getattr(time_merged_df,key+'_next'),
            t0,t1,
            t).values
        
    interp_snap['Coordinates'] = coords
    return interp_snap

def index_match_snapshots_with_dataframes(
    prev_sub_snap,
    next_sub_snap,
    extra_keys_to_extract=None):
    """
    keys_to_extract = ['Coordinates','Masses','SmoothingLength','ParticleIDs','ParticleChildIDsNumber']
    """
    
    init=time.time()
    print('Creating a merged DF')
    keys_to_extract = ['Coordinates','Masses','SmoothingLength','ParticleIDs','ParticleChildIDsNumber']
    if extra_keys_to_extract is not None:
        keys_to_extract += list(extra_keys_to_extract)

    ## convert snapshot dictionaries to pandas dataframes
    prev_df_snap = convertSnapToDF(prev_sub_snap,
        keys_to_extract=keys_to_extract)
    
    ## apparently index operations go faster if you sort by index
    prev_df_snap.sort_index(inplace=True)
    
    next_df_snap = convertSnapToDF(next_sub_snap,
        keys_to_extract=keys_to_extract)
    
    ## apparently index operations go faster if you sort by index
    next_df_snap.sort_index(inplace=True)
    
    
    ## remove particles that do not exist in the previous snapshot, 
    ##  difficult to determine which particle they split from
    next_df_snap_reduced = next_df_snap.reindex(prev_df_snap.index,copy=False)
    
    ## merge rows of dataframes based on 
    prev_next_merged_df_snap = pd.merge(
        prev_df_snap,
        next_df_snap_reduced,
        how='inner',
        on=prev_df_snap.index,
        suffixes=('','_next'),
        copy=False).set_index('key_0')
    
    ## remove particles that do not exist in the next snapshot, 
    ##  difficult to tell if they turned into stars or were merged. 
    prev_next_merged_df_snap = prev_next_merged_df_snap.dropna()
    
    return prev_next_merged_df_snap

def main():

    snapdir = "/scratch/projects/xsede/GalaxiesOnFIRE/metal_diffusion/m10q_res250/output"
    snapnum = 600 
    prev_galaxy = Galaxy(
        'm10q_res250',
        snapdir,
        snapnum,
        datadir='/scratch/04210/tg835099/data/metal_diffusion')

    times_gyr,snap_pairs = find_bordering_snapnums(
        prev_galaxy.snap_gyrs,
        dGyr=2)

    render_kwargs = {
        'weight_name':'Masses',
        'quantity_name':'Temperature',
        #min_quantity=2,
        #max_quantity=7,
        'quantity_adjustment_function':np.log10,
        'quick':True,
        'min_weight':-0.5,
        'max_weight':3,
        'weight_adjustment_function':lambda x: np.log10(x/(30**2/1200**2)) + 10 - 6, ## msun/pc^2,
        'cmap':'afmhot',
        'quick':True}

    control_flow(
        times_gyr,
        snap_pairs,
        render_kwargs=render_kwargs,
        name='m10q_res250',
        snapdir="/scratch/projects/xsede/GalaxiesOnFIRE/metal_diffusion/m10q_res250/output",
        datadir='/scratch/04210/tg835099/data/metal_diffusion')


if __name__ == '__main__':
    main()

