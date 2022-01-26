import numpy as np

from abg_python.plot_utils import plt
plt.rcParams['figure.dpi']=240
from abg_python.galaxy.gal_utils import ManyGalaxy

from firestudio.interpolate.interpolate import InterpolationHandler
from firestudio.studios.star_studio import StarStudio
from firestudio.studios.FIRE_studio import FIREStudio


def main(
    name='m12b_res7100',
    suite_name='metal_diffusion',
    multi_threads=1):


    many_galaxy = ManyGalaxy(name,suite_name=suite_name)

    last_galaxy = many_galaxy.loadAtSnapshot(many_galaxy.finsnap)
    snapnums,rcoms,rvirs = last_galaxy.get_rockstar_file_output()

    first_index = np.argmin(
        (many_galaxy.snapnums-snapnums[0])**2)

    interp_handler = InterpolationHandler(
        10, ## duration of movie, 10 sec
        many_galaxy.snap_gyrs[267], ## begininng time in Gyr
        many_galaxy.snap_gyrs[268]-0.001, ## end time in Gyr
        ## fixed position of camera, optionally could move camera around
        ## defines the fov as +- zdist
        camera_pos=[0,0,50],  
        ## how long should the line in the bottom left corner be?
        scale_line_length=10,
    )

    figs = interp_handler.interpolateAndRender(
        galaxy_kwargs={
            'ABG_force_multithread':8,
            ## flag to write snapshot w/ particles w/i rvir to disk
            'use_saved_subsnapshots':False, 
            'name':many_galaxy.name, ## 
            'final_orientation':True, ## face-on at z=0
            'loud_metadata':False, ## reduce print statements
            'use_rockstar_first':True, ## use rockstar over AHF if both exist
            'suite_name':suite_name},  ## path s.t. ~/snaps/{suite_name}/name/output
            ## use a soft-link:
            ## cd ~
            ## ln -s /scratch/projects/xsede/GalaxiesOnFIRE snaps
        render_kwargss=[
            {}, ## kwargs for FIREStudio render call
            {}], ## kwargs for StarStudio render call
        studio_kwargss=[
            {'savefig':None,'master_loud':True}, ## kwargs for FIREStudio initialization
            {'savefig':None,'maxden':2.2e8,'dynrange':4.7e2}], ## kwargs for StarStudio initialization
        multi_threads=multi_threads,
        which_studios=[FIREStudio,StarStudio],
        check_exists=True, ## skip rendering a frame if the png already exists
        timestamp=7.11419974, ## offset the timestamp by 0 Gyr ## bursty_time
        add_composition='og_test')  ## will add a composition frame of the requested Studios

if __name__ == '__main__':
    main()
