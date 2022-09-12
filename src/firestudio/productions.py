import warnings
import copy

import numpy as np

from .studios.gas_studio import GasStudio
from .studios.star_studio import StarStudio
from .studios.FIRE_studio import FIREStudio

class __Production(object):
    def __repr__(self):
        studio_str = repr(self.which_studio) if self.studio is None else repr(self.studio)

        return f"({self.name}) {studio_str}\n\tstudio_kwargs: {repr(self.studio_kwargs)}\n\trender_kwargs: {repr(self.render_kwargs)}"

    def __init__(self,name,production_kwargs): 
        self.name = name
        self.which_studio = production_kwargs['studio']
        self.render_kwargs = production_kwargs['render_kwargs']
        self.studio_kwargs = production_kwargs['studio_kwargs']
        self.studio = None
    
    def init(
        self,
        datadir:str,
        snapnum:int=None,
        name:str=None,
        **kwargs):

        for key,value in self.studio_kwargs.items():
            if key not in kwargs.keys(): kwargs[key] = value 

        self.studio = self.which_studio(datadir,snapnum,name,**kwargs)
    
    def render(self,ax,**kwargs):

        for key,value in self.render_kwargs.items():
            if key not in kwargs.keys(): 
                kwargs[key] = value if 'adjustment_function' not in key else lambda x: value(x,self.studio.Acell)
            


        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.studio.render(ax,**kwargs)

## not going to actually work, just needs to be a place-holder for now
def __msun_pc2_adjust_fn(x,Acell): 
    return np.log10(x/Acell) + 10 - 6 ## msun/pc^2

__mass_projection_kwargs = {
    'studio':GasStudio,
    'render_kwargs':{
        'weight_name':'Masses',
        'quantity_name':'Temperature',
        'cmap':'viridis',
        'min_weight':-0.5,
        'max_weight':3, 
        'weight_adjustment_function':__msun_pc2_adjust_fn
    },
    'studio_kwargs':{}
}

__velocity_projection_kwargs = {
    'studio':GasStudio,
    'render_kwargs':{
        'weight_name':'Masses',
        'quantity_name':'Vz',
        'cmap':'coolwarm',
        'min_quantity':-200,
        'max_quantity':200,
    },
    'studio_kwargs':{}
}

__temperature_projection_kwargs = {
    'studio':GasStudio,
    'render_kwargs':{
        'weight_name':'Masses',
        'quantity_name':'Temperature',
        'min_quantity':2,
        'max_quantity':7,
        'quantity_adjustment_function':np.log10,
    },
    'studio_kwargs':{}
}


__mock_hubble_kwargs = {
    'studio':StarStudio,
    'render_kwargs':{},
    'studio_kwargs':{'maxden':2.2e8,'dynrange':4.7e2} 
}

__fire_3_color_kwargs = {
    'studio':FIREStudio,
    'render_kwargs':{},
    'studio_kwargs':{} 
}

def register_production(prod_name,kwargs,loud=True):
    if loud: print('registering a production:',key[2:-len('_kwargs')]) 
    globals()[prod_name] = __Production(prod_name,kwargs)
    if loud: print(globals()[key[2:-len('_kwargs')]])

## algorithmically register 'productions' by the dictionaries that are listed above
for key in list(globals().keys()):
    if key[-len('_kwargs'):] == '_kwargs':
        register_production(key[2:-len('_kwargs')],globals()[key],loud=False)

## handle any productions which are perturbations on above productions
young_mock_hubble = copy.deepcopy(mock_hubble)
young_mock_hubble.render_kwargs.update({
    'no_dust':True,
    'age_max_gyr':25/1e3, ## 25 Myr
})