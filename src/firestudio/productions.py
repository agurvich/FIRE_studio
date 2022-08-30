import warnings

from .studios.gas_studio import GasStudio
from .studios.star_studio import StarStudio
from .studios.FIRE_studio import FIREStudio

class __Production(object):
    def __repr__(self):
        studio_str = repr(self.which_studio) if self.studio is None else repr(self.studio)

        return f"{studio_str}\n\tstudio_kwargs: {repr(self.studio_kwargs)}\n\trender_kwargs: {repr(self.render_kwargs)}"

    def __init__(self,production_kwargs): 
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
            if key not in kwargs.keys(): kwargs[key] = value 

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.studio.render(ax,**kwargs)

## not going to actually work, just needs to be a place-holder for now
def __msun_pc2_adjust_fn(x): return np.log10(x/my_studio.Acell) + 10 - 6 ## msun/pc^2

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
    globals()[prod_name] = __Production(kwargs)
    if loud: print(globals()[key[2:-len('_kwargs')]])

## algorithmically register 'productions' by the dictionaries that are listed above
for key in list(globals().keys()):
    if key[-len('_kwargs'):] == '_kwargs':
        register_production(key[2:-len('_kwargs')],globals()[key],loud=False)