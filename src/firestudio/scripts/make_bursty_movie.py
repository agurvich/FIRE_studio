import sys

import numpy as np

from abg_python.plot_utils import plt
plt.rcParams['figure.dpi']=240
from abg_python.galaxy.gal_utils import Galaxy

from firestudio.interpolate.interpolate import InterpolationHandler
from firestudio.studios.gas_studio import GasStudio
from firestudio.studios.star_studio import StarStudio
from firestudio.studios.composition import Composition


def main(
    coord_interp_mode='cylindrical',
    name='m12b_res57000',
    multi_threads=1):

    #galaxy = Galaxy('m12b_res7100',600)
    galaxy = Galaxy(name,600)
    _,bursty_time,_,_ = galaxy.get_bursty_regime(save_meta=True)

    interp_handler = InterpolationHandler(
        120,
        bursty_time-2,
        bursty_time+2,
        camera_pos=[0,0,30],
        scale_line_length=10,
        #time_slice=slice(-24*4,None),
        #time_slice=slice(-1,None),
        coord_interp_mode=coord_interp_mode,
    )

    figs = interp_handler.interpolateAndRender(
        {'name':galaxy.name,
        'keys_to_extract':['Metallicity','AgeGyr'],
        'final_orientation':True},
        render_kwargss=[
            {'weight_name':'Masses',
            'quantity_name':'Temperature',
            'min_quantity':2,
            'max_quantity':7,
            'quantity_adjustment_function':np.log10,
            #'min_weight':-0.5,
            #'max_weight':3, 
            },
            {'quick':False,#'use_metadata':True,
            'save_meta':True,#'assert_cached':False,
            #'min_quantity':2,'max_quantity':7,
            #'min_weight':-0.5,'max_weight':3 
            }], ## msun/pc^2,
        studio_kwargss=[
            {'figure_label':None},
            {'maxden':2.2e8,
            'dynrange':4.7e2,
            'no_dust':True,
            'scale_bar':False}],
        multi_threads=multi_threads,
        savefigs=None, ## defaults to GasStudio_frame and StarStudio_frame
        which_studios=[GasStudio,StarStudio],
        check_exists=False, ## skip rendering a frame if the png already exists
        timestamp=bursty_time,
        add_composition=True)  ## will add a composition frame of the requested Studios

if __name__ == '__main__':
    main()
