import copy
import numpy as np

from abg_python.function_utils import CLI_args
from abg_python.plot_utils import plt
plt.rcParams['figure.dpi']=300
from abg_python.galaxy.gal_utils import Galaxy

from firestudio.interpolate.interpolate import InterpolationHandler
from firestudio.studios.gas_studio import GasStudio
from firestudio.studios.star_studio import StarStudio
from firestudio.studios.FIRE_studio import FIREStudio
from firestudio.studios.composition import Composition

from firestudio.productions import mock_hubble,fire_3_color


@CLI_args()
def main(
    name='m12b_res7100',
    suite_name='metal_diffusion',
    multi_threads=10,
    take_avg_L=True,
    slice_index=None,
    nslices=None,
    duration=60, # s
    savefig='space_telescope'
    ):

    galaxy = Galaxy(name,600)
    _,bursty_time,_,_ = galaxy.get_bursty_regime(save_meta=True)
    galaxy.get_snapshotTimes()
    bursty_snap_index = np.argmin((galaxy.snap_gyrs-bursty_time)**2)

    if slice_index is not None and nslices is not None:
        nframes = duration * 24 ## default fps
        slice_width = nframes // nslices + ((nframes % nslices) != 0)
        time_slice = slice(slice_width*slice_index,slice_width*(slice_index+1))
        print(slice_index,time_slice,nframes)
    else: time_slice = None

    #time_slice = slice(0,24*6)

    interp_handler = InterpolationHandler(
        duration,
        galaxy.snap_gyrs[243]+1e-3,
        galaxy.snap_gyrs[396]-1e-3,
        camera_pos=[0,0,25],
        scale_line_length=5,
        snapshot_times=galaxy.snap_gyrs
        )

    ## select the images we want to make from predefined "productions"
    from firestudio.productions import mock_hubble,velocity_projection,fire_3_color,young_mock_hubble
    productions = [fire_3_color,mock_hubble]#velocity_projection,mock_hubble,young_mock_hubble]

    ## update the settings of the productions for our purposes here
    for production in productions:
        ## ensure we are doing a fresh projection every time (don't use the cache)
        ##  in this script but make it available for future use if we change our mind
        production.render_kwargs.update({'use_metadata':False,'save_meta':True})
        ## indicate that we want each image to get its own frame output (not just be included in the composition)
        production.studio_kwargs['savefig'] = f"{savefig}_{production.name}"

    figs = interp_handler.interpolateAndRender(
        {'name':galaxy.name,
        ##'ABG_force_multithread':60,
        'suite_name':suite_name,
        'take_avg_L':take_avg_L,
        'final_orientation':True,
        'loud_metadata':False},
         ## msun/pc^2,
        render_kwargss=[production.render_kwargs for production in productions], 
        studio_kwargss=[production.studio_kwargs for production in productions],
        multi_threads=multi_threads,
        which_studios=[production.which_studio for production in productions],
        check_exists=True, ## skip rendering a frame if the png already exists
        timestamp=bursty_time,
        add_composition=savefig,  ## will add a composition frame of the requested Studios
        time_slice=time_slice)

if __name__ == '__main__':
    main()
