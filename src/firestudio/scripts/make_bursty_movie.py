import sys
import getopt

import numpy as np

from abg_python.plot_utils import plt
plt.rcParams['figure.dpi']=300
from abg_python.galaxy.gal_utils import Galaxy

from firestudio.interpolate.interpolate import InterpolationHandler
from firestudio.studios.gas_studio import GasStudio
from firestudio.studios.star_studio import StarStudio
from firestudio.studios.composition import Composition


def main(
    name='m12b_res57000',
    multi_threads=1,
    take_avg_L = True 
    ):


    #galaxy = Galaxy('m12b_res7100',600)
    galaxy = Galaxy(name,600)
    _,bursty_time,_,_ = galaxy.get_bursty_regime(save_meta=True)
    galaxy.get_snapshotTimes()
    bursty_snap_index = np.argmin((galaxy.snap_gyrs-bursty_time)**2)


    interp_handler = InterpolationHandler(
        30,
        galaxy.snap_gyrs[bursty_snap_index]-1,
        galaxy.snap_gyrs[bursty_snap_index]-1e-3,#+2,
        camera_pos=[0,0,15],
        scale_line_length=5,
        #time_slice=slice(-24*4,None),
        #time_slice=slice(-10,None),
        time_slice=slice(0,100)
    )

    savefig = 'avg_L' if take_avg_L else 'rot_jhat'
    savefig+='_child'
    figs = interp_handler.interpolateAndRender(
        {'name':galaxy.name,
        ##'ABG_force_multithread':60,
        'take_avg_L':take_avg_L,
        'keys_to_extract':['Metallicity','AgeGyr'],
        'final_orientation':True,
        'loud_metadata':False},
        render_kwargss=[
            {'weight_name':'Masses',
            'quantity_name':'Temperature',
            'min_quantity':2,
            'max_quantity':7,
            'quantity_adjustment_function':np.log10,
            'save_meta':True,
            'use_metadata':False,
            #'assert_cached':False,
            #'min_weight':-0.5,
            #'max_weight':3, 
            },
            {'quick':False,#'use_metadata':True,
            'save_meta':True,#'assert_cached':False,
            'use_metadata':False,
            #'min_quantity':2,'max_quantity':7,
            #'min_weight':-0.5,'max_weight':3 
            },
            {
            #'save_meta':False,#'assert_cached':False,
            'use_metadata':False,
            }], ## msun/pc^2,
        studio_kwargss=[
            {
            'savefig':None
            },
            {'maxden':2.2e8,
            'dynrange':4.7e2,
            'savefig':None,
            #'no_dust':True,
            #'age_max_gyr':25/1e3, ## 25 Myr
            },
            {'savefig':'young_StarStudio',
            'maxden':2.2e8,
            'dynrange':4.7e2,
            'savefig':None,
            'no_dust':True,
            'age_max_gyr':25/1e3, ## 25 Myr
            }],
        multi_threads=multi_threads,
        which_studios=[GasStudio,StarStudio,StarStudio],
        check_exists=True, ## skip rendering a frame if the png already exists
        timestamp=bursty_time,
        add_composition=savefig)  ## will add a composition frame of the requested Studios

if __name__ == '__main__':
    argv = sys.argv[1:]
    opts,args = getopt.getopt(
        argv,'',[
        'name=',
        'multi_threads=',
        'take_avg_L='
        ])

    for i,opt in enumerate(opts):
        if opt[1]=='':
            opts[i]=('mode',opt[0].replace('-',''))
        else:
            try:
                ## if it's an int or a float this should work
                opts[i]=(opt[0].replace('-',''),eval(opt[1]))
            except:
                ## if it's a string... not so much
                opts[i]=(opt[0].replace('-',''),opt[1])
    main(**dict(opts))
