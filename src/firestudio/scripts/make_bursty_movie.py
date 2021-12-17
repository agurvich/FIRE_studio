import sys

import numpy as np

from abg_python.plot_utils import plt
plt.rcParams['figure.dpi']=240
from abg_python.galaxy.gal_utils import Galaxy

from firestudio.interpolate.interpolate import InterpolationHandler
from firestudio.studios.gas_studio import GasStudio


def main(coord_interp_mode='spherical'):
    galaxy = Galaxy('m12b_res7100',600)
    _,bursty_time,_,_ = galaxy.get_bursty_regime(save_meta=True)

    interp_handler = InterpolationHandler(
        120,
        bursty_time-2,
        bursty_time+2,
        camera_pos=[0,0,30],
        scale_line_length=10,
        #time_slice=slice(-24*8,None)
        coord_interp_mode=coord_interp_mode,
    )

    figs = interp_handler.interpolateAndRender(
        {'name':galaxy.name,
        'keys_to_extract':['Velocities','Temperature'],
        'final_orientation':True},
        render_kwargs={
            'quick':False,'use_metadata':False,
            'save_meta':True,
            'min_quantity':2,'max_quantity':7,
            'min_weight':-0.5,'max_weight':3 
            }, ## msun/pc^2,
    #    studio_kwargs={'maxden':2.2e8,'dynrange':4.7e2},
        multi_threads=7,
        savefig='two_color_%s'%coord_interp_mode,
        which_studio=GasStudio,
        check_exists=True) ## skip rendering a frame if the png already exists

if __name__ == '__main__':
    print(sys.argv)
    if len(sys.argv) != 2 or sys.argv[1] not in ['cylindrical','spherical']:
        raise KeyError("Must specify cylindrical or spherical as the only positional argument from the command line")
    main(sys.argv[1])
