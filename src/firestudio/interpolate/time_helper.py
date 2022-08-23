import gc
import itertools
import multiprocessing
import numpy as np
import pandas as pd
import os
import warnings

from abg_python.interpolate.time_interpolate_utils import convertToDF,cross_match_starformed_gas,finalize_df,make_interpolated_snap
from abg_python.galaxy.gal_utils import Galaxy
from abg_python.plot_utils import plt,ffmpeg_frames
from abg_python.physics_utils import addCircularVelocities

from firestudio.studios.gas_studio import GasStudio
from firestudio.studios.star_studio import StarStudio
from firestudio.studios.FIRE_studio import FIREStudio
from firestudio.studios.simple_studio import SimpleStudio
from firestudio.studios.composition import Composition

prev_galaxy,next_galaxy = None,None
merged_gas_df,merged_star_df = None,None

def single_threaded_control_flow(
    which_studios, ##Nstudios
    galaxy_kwargs,
    scene_kwargss, ## Ntimesteps
    studio_kwargss, ## Nstudios
    render_kwargss, ## Nstudios
    datadir,
    timestamp=0, ## offset by 0 Gyr
    polar=True):

    """ """

    if galaxy_kwargs is None: raise ValueError(
        'galaxy_kwargs must be a dictionary with,'+
        ' at minimum, name (i.e. m12i_res7100).')
    
    ## secret multithreading inside single threaded function, :eyes:
    if 'ABG_force_multithread' in galaxy_kwargs:
        ABG_force_multithread = galaxy_kwargs.pop('ABG_force_multithread')
    else: ABG_force_multithread = False
    
    ## determine what particle types we need to load from disk
    load_gas,load_star = get_load_flags(which_studios,render_kwargss)

    ## initialize the variables which we'll use
    ##  to avoid loading from disk unless we have to
    ##  i.e. moving next -> prev, etc...
    prev_snapnum,next_snapnum = None,None
    snapdict,star_snapdict = None,None

    snap_pairs = [scene_kwargs['snap_pair'] for scene_kwargs in scene_kwargss]
    if not ABG_force_multithread:
        return_values = []
        for scene_kwargs in scene_kwargss:     
            (prev_snapnum, next_snapnum,
            snapdict,star_snapdict,
            return_value) = render_this_scene(
                which_studios, ##Nstudios
                galaxy_kwargs, ## 1 dictionary shared by all scenes
                scene_kwargs, ## 1 dictionary for this scene
                studio_kwargss, ## Nstudios
                render_kwargss, ## Nstudios
                datadir,
                timestamp,
                load_gas,load_star,
                prev_snapnum,next_snapnum,
                snapdict,star_snapdict,
                ABG_force_multithread,
                polar)
            return_values +=[return_value]

    elif len(np.unique(snap_pairs)) == 2:
        pair = scene_kwargss[0]['snap_pair']
        (t0,t1) = scene_kwargss[0]['snap_pair_time']
        this_time = scene_kwargss[0]['time']
        snapdict,star_snapdict = get_interpolated_snaps(
            prev_snapnum,next_snapnum,
            pair,
            t0,t1,
            this_time,
            load_gas,load_star,
            polar=polar,
            **galaxy_kwargs)

        prev_snapnum,next_snapnum = pair
        ## set inside get_interpolated snaps, we won't need them
        ##  anymore
        global next_galaxy,prev_galaxy
        if hasattr(next_galaxy,'sub_snap'): del next_galaxy.sub_snap
        if hasattr(next_galaxy,'sub_star_snap'): del next_galaxy.sub_star_snap
        if hasattr(next_galaxy,'sub_star_snap'): del next_galaxy.sub_dark_snap
        if hasattr(prev_galaxy,'sub_snap'): del prev_galaxy.sub_snap
        if hasattr(prev_galaxy,'sub_star_snap'): del prev_galaxy.sub_star_snap
        if hasattr(prev_galaxy,'sub_star_snap'): del prev_galaxy.sub_dark_snap
        locals().keys()
        globals().keys()
        gc.collect()

        argss = zip(
            itertools.repeat(which_studios), ##Nstudios
            itertools.repeat(galaxy_kwargs), ## 1 dictionary shared by all scenes
            scene_kwargss, ## 1 dictionary for each scene
            itertools.repeat(studio_kwargss), ## Nstudios
            itertools.repeat(render_kwargss), ## Nstudios
            itertools.repeat(datadir),
            itertools.repeat(timestamp),
            itertools.repeat(load_gas),
            itertools.repeat(load_star),
            itertools.repeat(prev_snapnum),
            itertools.repeat(next_snapnum),
            itertools.repeat(snapdict),
            itertools.repeat(star_snapdict),
            itertools.repeat(ABG_force_multithread),
            itertools.repeat(polar))

        with multiprocessing.Pool(ABG_force_multithread) as pool:
            ## render_this_scene only returns FIRE studio return values
            ##  when ABG_force_multithread is true
            return_values = pool.starmap(render_this_scene,argss)
    else: raise ValueError("Not just a single pair of snapshots to use",np.unique(snap_pairs))

    return return_values

def render_this_scene(
    which_studios, ##Nstudios
    galaxy_kwargs, ## 1 dictionary shared by all scenes
    scene_kwargs, ## 1 dictionary for this scene
    studio_kwargss, ## Nstudios
    render_kwargss, ## Nstudios
    datadir,
    timestamp,
    load_gas,load_star,
    prev_snapnum,next_snapnum,
    snapdict,star_snapdict,
    ABG_force_multithread,
    polar=True):

    dummy_snap = {}
    dummy_snap['name'] = galaxy_kwargs['name']
    dummy_snap['datadir'] = datadir

    ## initialize dummy snapdict so that we can try and use
    ##  cached maps before loading anything from disk
    if 'snapnum' in galaxy_kwargs:
        dummy_snap['snapnum'] = galaxy_kwargs['snapnum']
        ## rely on user to provide figure_label explicitly
        timestamp = None 
    else:
        pair = scene_kwargs['snap_pair']
        dummy_snap['snapnum'] = pair[1]

    ## generate a timestamp if requested
    if timestamp is not None: 
        pair = scene_kwargs['snap_pair']
        (t0,t1) = scene_kwargs['snap_pair_time']
        scene_kwargs['figure_label'] = format_timestamp(
            pair,
            scene_kwargs['time'],
            t0,t1,
            timestamp)

    ## attempt to use cached images
    try:
        ## call the abstract draw function
        return_value = multi_worker_function(
            which_studios,
            dummy_snap,
            dummy_snap,
            scene_kwargs,
            studio_kwargss,
            [{'assert_cached':True,'loud':False,**render_kwargs} for render_kwargs in render_kwargss])
        if not ABG_force_multithread:
            return prev_snapnum,next_snapnum,snapdict,star_snapdict,return_value
        else: return return_value
    ## yeah alright, it was worth a shot. let's load from disk and project then
    except (AssertionError,KeyError): pass #raise

    ## determine if we've already loaded the data or if we need to open it from disk
    snapdict,star_snapdict,prev_snapnum,next_snapnum = load_data_from_disk_if_necessary(
        galaxy_kwargs,
        scene_kwargs,
        prev_snapnum,next_snapnum,
        snapdict,star_snapdict,
        load_gas,load_star,
        polar=polar) ## determines what coordinates we save in snapdict

    ## call the abstract draw function with live data this time
    return_value = multi_worker_function(
        which_studios,
        snapdict,
        star_snapdict,
        scene_kwargs,
        studio_kwargss,
        render_kwargss)

    if not ABG_force_multithread:
        return prev_snapnum,next_snapnum,snapdict,star_snapdict,return_value
    else: return return_value

def render_ffmpeg_frames(studio_kwargss,galaxy_kwargs,nframes,fps):
    for studio_kwargs in studio_kwargss:
        if 'keys_to_extract' in galaxy_kwargs: galaxy_kwargs.pop('keys_to_extract')
        if 'use_saved_subsnapshots' in galaxy_kwargs: galaxy_kwargs.pop('use_saved_subsnapshots')
        if 'force_theta_TB' in galaxy_kwargs: galaxy_kwargs.pop('force_theta_TB')
        if 'force_phi_TB' in galaxy_kwargs: galaxy_kwargs.pop('force_phi_TB')

        if 'snapnum' not in galaxy_kwargs: galaxy_kwargs['snapnum'] = None
        galaxy = Galaxy(**galaxy_kwargs)

        this_savefig = studio_kwargs['savefig']
        ## in case we rendered multiple types of frames
        ##  we'll loop through a list of savefigs (even if there's only one)
        if this_savefig is not None:
            format_str = '%s'%this_savefig + '_frame_%0'+'%dd.png'%(np.ceil(np.log10(nframes)))
            ## ffmpeg the frames
            ffmpeg_frames(
                os.path.join(galaxy.datadir,'firestudio'),
                [format_str],
                savename=galaxy_kwargs['name'],
                framerate=fps,
                extension='.mp4')

def load_data_from_disk_if_necessary(
    galaxy_kwargs,
    scene_kwargs,
    prev_snapnum,next_snapnum,
    snapdict,star_snapdict,
    load_gas,load_star,
    polar=True):

    ## determine if we are using 'just' a single snapshot's data
    if 'snapnum' in galaxy_kwargs: 
        if snapdict is None:
            ## only open the snapshot data once
            snapdict,star_snapdict = get_single_snap(
                load_gas,load_star,
                **galaxy_kwargs)
    ## or if we want to interpolate between two snapshots
    else:
        ## remove timing information for this dictionary
        pair = scene_kwargs['snap_pair']
        (t0,t1) = scene_kwargs['snap_pair_time']
        this_time = scene_kwargs['time']
        snapdict,star_snapdict = get_interpolated_snaps(
            prev_snapnum,next_snapnum,
            pair,
            t0,t1,
            this_time,
            load_gas,load_star,
            polar=polar,
            **galaxy_kwargs)
        prev_snapnum,next_snapnum = pair
    return snapdict,star_snapdict,prev_snapnum,next_snapnum

def get_single_snap(
    load_gas,
    load_star,
    force_theta_TB=None,
    force_phi_TB=None,
    keys_to_extract=None,
    use_saved_subsnapshots=True,
    extract_DM=True, ## ignored
    **galaxy_kwargs):

    if keys_to_extract is None: keys_to_extract = []

    snapdict,star_snapdict = {},{}

    galaxy = Galaxy(**galaxy_kwargs)
    if not extract_DM: 
        print('ignoring user request to not extract DM b.c.'
            +' we need it for Vc which is used in the interpolation.')
        extract_DM=True
    galaxy.extractMainHalo(
        compute_stellar_hsml=load_star,
        force_theta_TB=force_theta_TB,
        force_phi_TB=force_phi_TB,
        use_saved_subsnapshots=use_saved_subsnapshots,
        save_meta=use_saved_subsnapshots,
        loud=False,
        extract_DM=extract_DM, ## force to be true
        jhat_coords=False)

    if load_gas: 
        for key in galaxy.sub_snap.keys():
            if key in keys_to_extract:
                snapdict[key] = galaxy.sub_snap[key] 
        snapdict['name'] = galaxy.name
        snapdict['datadir'] = galaxy.datadir
        snapdict['snapnum'] = galaxy.snapnum

    if load_star: 
        star_snapdict = galaxy.sub_star_snap
        for key in galaxy.sub_star_snap.keys():
            if key in keys_to_extract:
                star_snapdict[key] = galaxy.sub_star_snap[key] 
        star_snapdict['name'] = galaxy.name
        star_snapdict['datadir'] = galaxy.datadir
        star_snapdict['snapnum'] = galaxy.snapnum

    ## add the circular velocity for each particle by computing
    ##  the potential energy it has using a look-up table
    addCircularVelocities(
        [galaxy.sub_snap,galaxy.sub_star_snap,galaxy.sub_dark_snap],
        [snapdict]*load_gas + [star_snapdict]*load_star)


    del galaxy

    return snapdict,star_snapdict
    
def get_interpolated_snaps(
    prev_snapnum,
    next_snapnum,
    pair,
    t0,t1,
    this_time,
    load_gas,load_star,
    keys_to_extract=None,
    polar=True, ## use polar interpolation for coordinates?
    take_avg_L=False,
    **galaxy_kwargs):

    if keys_to_extract is None: keys_to_extract = []

    global prev_galaxy,next_galaxy
    global merged_gas_df,merged_star_df
    #if not i%10: print(i,pair,pair_times)
    ## determine if the galaxies in the pair are actually
    ##  changed, and if so, open their data from the disk.
    prev_galaxy,next_galaxy,changed = load_gals_from_disk(
        prev_snapnum,next_snapnum,
        pair,
        prev_galaxy,
        next_galaxy,
        compute_stellar_hsml=load_star,
        polar=polar,
        **galaxy_kwargs)

    ## update the previous/next snapnums
    prev_snapnum,next_snapnum = pair
    ## make an interpolated snapshot with these galaxies,
    ##  this takes a while so we'll hold onto it and only 
    ##  make a new one if necessary.
    #changed = True
    if changed:
        ## create the gas dataframes if necessary
        if load_gas: 
            #print("converting gas to DF")
            prev_gas_df = convertToDF(prev_galaxy.sub_snap,keys_to_extract,polar)
            next_gas_df = convertToDF(next_galaxy.sub_snap,keys_to_extract,polar)
        
        ## create the star dataframes if necessary
        if load_star:
            #print("converting stars to DF")
            prev_star_df = convertToDF(prev_galaxy.sub_star_snap,keys_to_extract,polar)
            next_star_df = convertToDF(next_galaxy.sub_star_snap,keys_to_extract,polar)

        ## if we have both, we can try to cross-match the positions
        ##  of particles between them
        if load_gas and load_star:
            #print("cross-matching starformed gas")

            ## add AgeGyr arrays which will be used to disappear gas particles
            ##  when we multi index match them with splits/stars
            prev_gas_df['AgeGyr'] = np.ones(prev_gas_df.shape[0])
            next_gas_df['AgeGyr'] = np.ones(next_gas_df.shape[0])

            (prev_gas_df,
            next_gas_df,
            prev_star_df,
            next_star_df) = cross_match_starformed_gas(
                t0,t1,
                prev_gas_df,next_gas_df,
                prev_star_df,next_star_df)


        if load_gas:
            ## look for ancestor gas particles for those particles which split
            """
            merged_gas_df = search_multi_ids(
                merged_gas_df, ## df to search for ancestors
                merged_gas_df, ## df w/ targets
                forward=False) ## look in the past
                
            if load_star:
                ## find final position that gas particles which turn into stars
                ##  should interpolate toward
                search_multi_ids(
                    merged_star_df, ## df to use as lookup
                    merged_gas_df, ## df w/ targets
                    forward=True) ## look in the future
            """
            ## merge rows of dataframes based on particle ID
            merged_gas_df = prev_gas_df.join(
                next_gas_df,
                rsuffix='_next',
                how='outer')
            
            del prev_galaxy.sub_snap
            #del next_galaxy.sub_snap ## don't delete this b.c. we'll need it when we load the next one

            ## fill values w. extrapolation in both directions
            ##  add polar coordinates and velocities, and drop any remaining nans or duplicates
            merged_gas_df = finalize_df(
                t0,t1,
                merged_gas_df,
                polar=polar,
                take_avg_L=take_avg_L)

        if load_star: 
            ## merge rows of dataframes based on particle ID
            merged_star_df = prev_star_df.join(
                next_star_df,
                rsuffix='_next',
                how='outer')

            del prev_galaxy.sub_star_snap
            #del next_galaxy.sub_star_snap ## don't delete this b.c. we'll need it when we load the next one

            ## fill values w. extrapolation in both directions
            ##  add polar coordinates and velocities, and drop any remaining nans or duplicates
            merged_star_df = finalize_df(
                t0,t1,
                merged_star_df,
                polar=polar,
                take_avg_L=take_avg_L)
    else: pass ## [ if changed: ] 

    ## create the interp_snap with new values for the new time
    if load_gas: interp_gas_snapdict = make_interpolated_snap(
        t0,t1,
        merged_gas_df,
        this_time) 
    else: interp_gas_snapdict = {}            
    if load_star: interp_star_snapdict = make_interpolated_snap(
        t0,t1,
        merged_star_df,
        this_time)
    else: interp_star_snapdict = {}

    ## keep outside the conditional b.c. worker function
    ##  looks for them in gas_snapdict
    interp_gas_snapdict['name'] = next_galaxy.name
    interp_gas_snapdict['datadir'] = next_galaxy.datadir
    interp_gas_snapdict['snapnum'] = next_galaxy.snapnum
    interp_gas_snapdict['this_time_Gyr'] = this_time
    interp_gas_snapdict['prev_time_Gyr'] = t0
    interp_gas_snapdict['next_time_Gyr'] = t1

    interp_star_snapdict['name'] = next_galaxy.name
    interp_star_snapdict['datadir'] = next_galaxy.datadir
    interp_star_snapdict['snapnum'] = next_galaxy.snapnum
    interp_star_snapdict['this_time_Gyr'] = this_time
    interp_star_snapdict['prev_time_Gyr'] = t0
    interp_star_snapdict['next_time_Gyr'] = t1


    return interp_gas_snapdict,interp_star_snapdict


def load_gals_from_disk(
    prev_snapnum,next_snapnum,
    pair,
    prev_galaxy,next_galaxy,
    testing=False,
    compute_stellar_hsml=False,
    force_theta_TB=None,
    force_phi_TB=None,
    use_saved_subsnapshots=True,
    extract_DM=True,
    polar=True,
    **kwargs):
    """ Determines whether it needs to load a new galaxy from disk
        or if we already have what we need."""

    if not extract_DM: 
        print('ignoring user request to not extract DM b.c.'
            +' we need it for Vc which is used in the interpolation.')
        extract_DM=True

    ## -- check the prev galaxy
    ## keep the current snapnum
    if pair[0] == prev_snapnum: prev_galaxy=prev_galaxy
    ## step forward in time, swap pointers
    elif pair[0] == next_snapnum:
        prev_galaxy = next_galaxy
        next_galaxy = None
    ## will need to read from disk
    else: prev_galaxy = None
    
    ## -- now the next galaxy
    ## keep the current snapnum
    if pair[1] == next_snapnum: next_galaxy = next_galaxy
    ## will need to read from disk
    else: next_galaxy = None

    changed = False ## flag for whether we loaded something from disk
    if prev_galaxy is None:
        print('loading',pair[0],'from disk')
        if not testing:
            prev_galaxy = Galaxy(snapnum=pair[0],**kwargs)
            prev_galaxy.extractMainHalo(
                compute_stellar_hsml=compute_stellar_hsml,
                force_theta_TB=force_theta_TB,
                force_phi_TB=force_phi_TB,
                use_saved_subsnapshots=use_saved_subsnapshots,
                save_meta=use_saved_subsnapshots,
                loud=False,
                extract_DM=extract_DM,
                jhat_coords=False)
            if polar: 
                addCircularVelocities(
                    [prev_galaxy.sub_snap,
                    prev_galaxy.sub_star_snap,
                    prev_galaxy.sub_dark_snap],
                    [prev_galaxy.sub_snap,
                    prev_galaxy.sub_star_snap])
        else: prev_galaxy = pair[0]
        changed = True
    if next_galaxy is None:
        print('loading',pair[1],'from disk')
        if not testing:
            next_galaxy = Galaxy(snapnum=pair[1],**kwargs)
            next_galaxy.extractMainHalo(
                compute_stellar_hsml=compute_stellar_hsml,
                force_theta_TB=force_theta_TB,
                force_phi_TB=force_phi_TB,
                use_saved_subsnapshots=use_saved_subsnapshots,
                save_meta=use_saved_subsnapshots,
                loud=False,
                extract_DM=extract_DM,
                jhat_coords=False)
            if polar: 
                addCircularVelocities(
                    [next_galaxy.sub_snap,
                    next_galaxy.sub_star_snap,
                    next_galaxy.sub_dark_snap],
                    [next_galaxy.sub_snap,
                    next_galaxy.sub_star_snap])
        else: next_galaxy = pair[1]
        changed = True
        
    return prev_galaxy,next_galaxy,changed

def format_timestamp(snap_pair,t,t0,t1,offset=0):
    ## Myr precision
    this_string = "%d - %.3f Gyr - %d Myr - %d Myr"%(snap_pair[0],t-offset,(t-t0)*1e3,(t1-t)*1e3)

    return this_string

def get_load_flags(which_studios,render_kwargss):
    ## determine which data we'll have to open
    load_gas,load_star = False,False
    for which_studio,render_kwargs in zip(which_studios,render_kwargss):
        this_load_gas,this_load_star = load_data_flags(which_studio,render_kwargs)
        load_gas = load_gas or this_load_gas
        load_star = load_star or this_load_star
    return load_gas,load_star

def load_data_flags(which_studio,render_kwargs):
    ## determine which data needs to be loaded into shared memory
    if which_studio is not GasStudio:
        #which_studio is StarStudio or which_studio is FIREStudio:
        load_gas = True
        load_star = True
    elif 'snapdict_name' in render_kwargs and render_kwargs['snapdict_name'] == 'star':
        load_gas = False
        load_star = True
    else:
        load_gas = True
        load_star = False
    return load_gas,load_star

def multi_worker_function(
    which_studios,
    this_snapdict,
    this_star_snapdict=None,
    scene_kwargs=None,
    studio_kwargss=None,
    add_render_kwargss=None):

    if len(which_studios) != len(studio_kwargss) != len(add_render_kwargss):
        raise ValueError(
            "Must have matching length lists of: "+
            "studios, studio kwargs, and render kwargs")

    for which_studio,studio_kwargs,add_render_kwargs in zip(
        which_studios,studio_kwargss,add_render_kwargss):
        try:
            worker_function(
                which_studio,
                this_snapdict,
                this_star_snapdict,
                {**scene_kwargs,**studio_kwargs},
                add_render_kwargs)
        except AssertionError: raise
        except Exception as e:
            print(
                '%03d'%(this_snapdict['snapnum']-1),
                '%03d'%this_snapdict['snapnum'],
                'failed',
                e.__context__)
            raise

def worker_function(
    which_studio,
    this_snapdict,
    this_star_snapdict=None,
    studio_kwargs=None,
    add_render_kwargs=None):

    if studio_kwargs is None: studio_kwargs = {'savefig':None}
    if add_render_kwargs is None: add_render_kwargs = {}

    if studio_kwargs['savefig'] is not None and 'savefig_suffix' in studio_kwargs:
        this_savefig = studio_kwargs['savefig'] + studio_kwargs.pop('savefig_suffix')
    else: this_savefig = None

    ## decide what we want to pass to the GasStudio
    if which_studio is GasStudio: render_kwargs = {
        'weight_name':'Masses',
        'quantity_name':'Temperature',
        'min_quantity':2,
        'max_quantity':7,
        'quantity_adjustment_function':np.log10,
        #'min_weight':-0.5,
        #'max_weight':3,
        }
    elif which_studio is StarStudio: render_kwargs = {}
    elif which_studio is FIREStudio: render_kwargs = {}
    elif which_studio is SimpleStudio: render_kwargs = {}
    elif which_studio is Composition: 
        studios_tuple = studio_kwargs.pop('studios_tuple')
        which_studio_actual = which_studio
        which_studio = lambda *args,**kwargs: which_studio_actual(
            studios_tuple,*args,**kwargs)
        render_kwargs = {}
    else: raise TypeError("%s is not GasStudio or StarStudio"%repr(which_studio))

    render_kwargs.update(add_render_kwargs)

    #print(which_studio)
    #print(studio_kwargs)
    #print(render_kwargs)
    if ('time' in studio_kwargs.keys() and studio_kwargs['time'] is not None): 
        setup_id_append = "_time%.5f"%studio_kwargs['time']
    else: setup_id_append = ''

    my_studio = which_studio(
        os.path.join(this_snapdict['datadir'],'firestudio'),
        this_snapdict['snapnum'], ## attribute this data to the next_snapnum's projection file
        this_snapdict['name'],
        gas_snapdict=this_snapdict,
        star_snapdict=this_star_snapdict,
        setup_id_append=setup_id_append,
        **{'master_loud':False,**studio_kwargs,'savefig':this_savefig})

    if which_studio is GasStudio and render_kwargs['weight_name'] == 'Masses':
        render_kwargs['weight_adjustment_function'] = lambda x: np.log10(x/my_studio.Acell) + 10 - 6 ## msun/pc^2,
    

    ## ignore warnings to reduce console spam
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ## overwrite loud rendering to reduce console spam
        ax,im = my_studio.render(None,**{**render_kwargs,'loud':False})
    axs = np.array([ax]).reshape(-1)
    fig = axs[0].get_figure()

    if my_studio.savefig is not None: 
        plt.close(fig)
        print(
            '%03d'%(this_snapdict['snapnum']-1),
            '%03d'%this_snapdict['snapnum'],
            my_studio.savefig)
    else: return fig
