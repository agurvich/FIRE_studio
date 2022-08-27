import numpy as np
import sys

from abg_python.plot_utils import plt
plt.rcParams['figure.dpi']=300
from abg_python.galaxy.gal_utils import ManyGalaxy

from firestudio.interpolate.interpolate import InterpolationHandler
from firestudio.studios.star_studio import StarStudio
from firestudio.studios.FIRE_studio import FIREStudio
from firestudio.studios.simple_studio import SimpleStudio 


def main(
    name='m12b_res7100',
    suite_name='metal_diffusion',
    multi_threads=4,
    savefig_str=''):

    if savefig_str == '': savefig_str = 'multi_cam_test'


    #many_galaxy = ManyGalaxy(name,suite_name=suite_name)
    #last_galaxy = many_galaxy.loadAtSnapshot(many_galaxy.finsnap)
    #last_galaxy.get_snapshotTimes()
    #snapnums,rcoms,rvirs = last_galaxy.get_rockstar_file_output()

    #first_index = np.argmin((last_galaxy.snapnums-snapnums[0])**2)

    interp_handler = InterpolationHandler(
        30, ## duration of movie, 10 sec
        ## fixed position of camera, optionally could move camera around
        ## defines the fov as +- zdist
        camera_pos=[0,0,50],  
        ## how long should the line in the bottom left corner be?
        scale_line_length=10,
    )

    interp_handler.scene_handler.add_keyframe(
        interp_handler.scene_handler.total_duration_sec, ## just have 2 keyframes, put this one at the end
        camera_pos=[0,0,25],  
        scale_line_length=5)

    figs = interp_handler.interpolateAndRender(
        galaxy_kwargs={
            #'ABG_force_multithread':10,
            ## flag to write snapshot w/ particles w/i rvir to disk
            'snapnum':600,
            #'take_avg_L':True,
            'use_saved_subsnapshots':True, 
            'name':name, ## 
            'final_orientation':True, ## face-on at z=0
            'loud_metadata':False, ## reduce print statements
            'suite_name':suite_name},  ## path s.t. ~/snaps/{suite_name}/name/output
            ## use a soft-link:
            ## cd ~
            ## ln -s /scratch/projects/xsede/GalaxiesOnFIRE snaps
        render_kwargss=[
            {
                #'age_max_gyr':25/1e3, ## 25 Myr
                'use_metadata':False,
                'save_meta':False
            }], ## kwargs for StarStudio render call
        studio_kwargss=[
            {
                'savefig':savefig_str,
                # 'maxden':2.2e8,'dynrange':4.7e2,
            #'no_dust':True,
            }], ## kwargs for StarStudio initialization
        multi_threads=None,
        shared_memory=True,
        which_studios=[FIREStudio],
        check_exists=True, ## skip rendering a frame if the png already exists
        #timestamp=last_galaxy.get_bursty_regime()[0]/1e3, ## offset the timestamp by 0 Gyr ## 
        add_composition=False)  ## will add a composition frame of the requested Studios

if __name__ == '__main__':
    main(savefig_str=sys.argv[1])
