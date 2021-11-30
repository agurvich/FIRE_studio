import numpy as np

from abg_python.interpolate.time_interpolate_utils import index_match_snapshots_with_dataframes,make_interpolated_snap
from abg_python.galaxy.gal_utils import Galaxy

#### MPS LOAD BALANCING
def split_into_n_approx_equal_chunks(snap_pairs,nchunks):
    
    nrenders = len(snap_pairs)
    ## split matching pairs into groups
    indices = split_pairs(snap_pairs[:])
    splitted = np.array_split(snap_pairs,indices)
    
    ## determine how many of each matching pair there are
    n_renders_per_split = [len(this_split) for this_split in splitted]
    
    mps_chunks = []
    this_chunk = 0
    # 1 would bias early, do this if the remainder would dump a lot of pairs on the last
    ##  process, i.e. when the remainder > nchunks/2
    per_chunk = nrenders//nchunks + (nrenders%nchunks > nchunks//2)
    
    for i in range(len(indices)+1):
        #print(this_chunk,per_chunk)
        if this_chunk >= per_chunk:
            mps_chunks+=[indices[i-1]]
            this_chunk=0
        this_chunk+=n_renders_per_split[i]
    
    mps_chunks = np.array(mps_chunks)#-indices[0]
    mps_chunks=list(mps_chunks)
    print('split into:',np.diff([0]+mps_chunks+[len(snap_pairs)]))   
    return mps_chunks

def split_pairs(snap_pairs):
    indices = []
    changed = split_head(snap_pairs)
    cur_index = 0
    while changed > 0:
        changed = split_head(snap_pairs[cur_index:])
        cur_index+=changed
        indices+=[cur_index]

    return indices[:-1]#np.array_split(snap_pairs,indices[:-1])
    
def split_head(snap_pairs):
    split_index = np.argmax(np.logical_not(np.all(snap_pairs==snap_pairs[0],axis=1)))
    return split_index

def single_threaded_control_flow(
    which_studio,
    times,
    snap_pairs,
    snap_pair_times,
    galaxy_kwargs,
    studio_kwargs,
    render_kwargs,
    savefigs):
    """ """

    ## put here to avoid circular import
    from .interpolate import worker_function

    if galaxy_kwargs is None: raise ValueError(
        'galaxy_kwargs must be a dictionary with,'+
        ' at minimum, name (i.e. m12i_res7100).')

    prev_galaxy,next_galaxy = None,None
    prev_snapnum,next_snapnum = None,None

    return_values = []

    extra_keys_to_extract = (
        galaxy_kwargs.pop('extra_keys_to_extract') if 
        'extra_keys_to_extract' in galaxy_kwargs else [])

    for i,(pair,pair_times) in enumerate(zip(snap_pairs,snap_pair_times)):     

        savefig = savefigs[i]

        if not i%10: print(i,pair,pair_times)
        ## determine if the galaxies in the pair are actually
        ##  changed, and if so, open their data from the disk.
        prev_galaxy,next_galaxy,changed = load_gals_from_disk(
            prev_snapnum,next_snapnum,
            pair,
            prev_galaxy,
            next_galaxy,
            **galaxy_kwargs)

        ## update the previous/next snapnums
        prev_snapnum,next_snapnum = pair

        this_time = times[i]
        
        if changed:
            ## make an interpolated snapshot with these galaxies,
            ##  this takes a while so we'll hold onto it and only 
            ##  make a new one if necessary.
            t0,t1 = pair_times
            time_merged_df = index_match_snapshots_with_dataframes(
                prev_galaxy.sub_snap,
                next_galaxy.sub_snap,
                extra_keys_to_extract=extra_keys_to_extract)

        ## update the interp_snap with new values for the new time
        interp_snap = make_interpolated_snap(this_time,time_merged_df,t0,t1)
        ## TODO
        print(interp_snap.keys(),"TODO in time_helper, figure out what keys are in interp_snap")

        ## TODO interpolate on stars as well
        interp_star_snap = None

        ## call the function we were passed
        return_values += [worker_function(
            which_studio,
            interp_snap,
            interp_star_snap,
            studio_kwargs,
            render_kwargs,
            savefig)]

    return return_values

def load_gals_from_disk(
    prev_snapnum,next_snapnum,
    pair,
    prev_galaxy,next_galaxy,
    testing=False,
    **kwargs):
    """ Determines whether it needs to load a new galaxy from disk
        or if we already have what we need."""

    ## -- check the prev galaxy
    ## keep the current snapnum
    if pair[0] == prev_snapnum:
        prev_galaxy=prev_galaxy
    ## step forward in time, swap pointers
    elif pair[0] == next_snapnum:
        prev_galaxy = next_galaxy
        next_galaxy = None
    ## will need to read from disk
    else:
        prev_galaxy = None
    
    ## -- now the next galaxy
    ## keep the current snapnum
    if pair[1] == next_snapnum:
        next_galaxy = next_galaxy
    ## will need to read from disk
    else:
        next_galaxy = None

    changed = False ## flag for whether we loaded something from disk
    if prev_galaxy is None:
        print('loading',pair[0],'from disk')
        if not testing:
            prev_galaxy = Galaxy(snapnum=pair[0],**kwargs)
            prev_galaxy.extractMainHalo()
        else: prev_galaxy = pair[0]
        changed = True
    if next_galaxy is None:
        print('loading',pair[1],'from disk')
        if not testing:
            next_galaxy = Galaxy(snapnum=pair[1],**kwargs)
            next_galaxy.extractMainHalo()
        else: next_galaxy = pair[1]
        changed = True
        
    return prev_galaxy,next_galaxy,changed