import numpy as np
import sys

from abg_python.plot_utils import plt
plt.rcParams['figure.dpi']=300
from abg_python.galaxy.gal_utils import ManyGalaxy

from firestudio.interpolate.interpolate import InterpolationHandler,SceneInterpolationHandler


def set_random_orientations(
    scene_handler:SceneInterpolationHandler,
    radius,
    nframes=None):

    scene_handler.add_keyframe(
        None,
        camera_pos=[radius,0,0], ## add edge-on view
        nsteps=1,
        loud=False)

    if nframes is None: nframes = scene_handler.nframes - len(scene_handler.scene_kwargss)
    ## generate random orientations
    np.random.seed(51694) ## it's my birthday :3
    phis = np.random.uniform(0,2*np.pi,size=nframes)
    costhetas = np.random.uniform(-1,1,size=nframes)
    thetas = np.arccos(costhetas) ## goes from 0 -> pi, just like we want
    camera_poss = np.array([
        np.sin(thetas)*np.cos(phis),
        np.sin(thetas)*np.sin(phis),
        costhetas]).T
    
    radii = np.random.uniform(0.5,1,size=nframes)*radius

    for radius,camera_pos in zip(radii,camera_poss):
        scene_handler.add_keyframe(
            None,
            camera_pos=radius*camera_pos,  
            nsteps=1,
            loud=False)

def main(
    name='m12b_res7100',
    suite_name='metal_diffusion',
    savefig_str='',
    radius=50,
    pixels=25,
    snapnum=243):

    if savefig_str == '': savefig_str = 'multi_cam_test'


    many_galaxy = ManyGalaxy(name,suite_name=suite_name)
    last_galaxy = many_galaxy.loadAtSnapshot(many_galaxy.finsnap)
    last_galaxy.get_snapshotTimes()
    #snapnums,rcoms,rvirs = last_galaxy.get_rockstar_file_output()

    #first_index = np.argmin((last_galaxy.snapnums-snapnums[0])**2)

    interp_handler = InterpolationHandler(
        1, ## duration of movie, 1 sec
        ## fixed position of camera, optionally could move camera around
        ## defines the fov as +- zdist
        camera_pos=[0,0,radius],   ## face-on
        ## how long should the line in the bottom left corner be?
        scale_line_length=10,
        fps=24 ## 1 sec + 24 fps -> 24 orientations
    )
 
    set_random_orientations(
        interp_handler.scene_handler,
        radius)

    from firestudio.productions import mock_hubble,velocity_projection

    mock_hubble.studio_kwargs.update({
        'figure_label':f'{snapnum} - redshift = {last_galaxy.snap_zs[snapnum]:0.2f}',
        'savefig':savefig_str+'_star',
        'pixels':pixels,
    })
    mock_hubble.render_kwargs.update({
        #'no_dust':True,
        #'age_max_gyr':25/1e3, ## 25 Myr
        #'use_metadata':False,
        'save_meta':True,
        'inset_axes':True
    })

    velocity_projection.studio_kwargs.update({
        'savefig':savefig_str+'_gas',
        'pixels':pixels,
    })
    velocity_projection.render_kwargs.update({
        #'use_metadata':False,
        'save_meta':True
    })


    productions = [mock_hubble,velocity_projection]
    figs = interp_handler.interpolateAndRender(
        galaxy_kwargs={
            #'ABG_force_multithread':10,
            ## flag to write snapshot w/ particles w/i rvir to disk
            'snapnum':snapnum,
            #'take_avg_L':True,
            'use_saved_subsnapshots':True, 
            'name':name, ## 
            'final_orientation':True, ## face-on at z=0
            'loud_metadata':False, ## reduce print statements
            'suite_name':suite_name},  ## path s.t. ~/snaps/{suite_name}/name/output
            ## use a soft-link:
            ## cd ~
            ## ln -s /scratch/projects/xsede/GalaxiesOnFIRE snaps
        render_kwargss=[production.render_kwargs for production in productions], 
        studio_kwargss=[production.studio_kwargs for production in productions],
        multi_threads=None,
        shared_memory=True,
        which_studios=[production.which_studio for production in productions],
        check_exists=True, ## skip rendering a frame if the png already exists
        #timestamp=last_galaxy.get_bursty_regime()[0]/1e3, ## offset the timestamp by 0 Gyr ## 
        add_composition=False)  ## will add a composition frame of the requested Studios

if __name__ == '__main__':
    main(savefig_str=sys.argv[1])
