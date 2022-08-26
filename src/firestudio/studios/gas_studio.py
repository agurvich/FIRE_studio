## builtin imports
import os
import numpy as np 
import ctypes
import copy


## abg_python imports
from abg_python import append_function_docstring,append_string_docstring
from abg_python.plot_utils import plt,addColorbar,get_cmap
from abg_python.galaxy.metadata_utils import metadata_cache

## firestudio imports
from .studio import Studio

def vfloat(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

class GasStudio(Studio):
    """ `FIREstudio` class for making gas projection images.
        Can be used for stars, but you will either have to pass smoothing lengths
        in or allow FIREstudio to calculate them itself, which  can take a long time. 

Important methods include: 

* [`GasStudio.weightAvgAlongLOS`](#gasstudioweightavgalonglos) 
* [`GasStudio.render`](#gasstudiorender) 
* [`Studio.__init__`](#studio__init__) 
* [`Studio.set_ImageParams`](#studioset_imageparams)"""


    required_snapdict_keys = ['Coordinates','Masses','SmoothingLength','Temperature','Velocities']

    def __repr__(self):
        """ implementation of built-in __repr__ method for printing.

        Returns
        -------
        str
            string to print back to the console
        """
        return 'GasStudio instance'

    def set_ImageParams(
        self,
        use_defaults:bool=False,
        loud:bool=True,
        **kwargs):
        """ Changes the parameters of the image such as camera orientation, frame size, etc. 
            If `use_defaults=True` then default values of the parameters will be set and will 
            overwrite the current state. Leave `use_defaults=False` to adjust only the keywords passed. 

        Parameters
        ----------
        use_defaults : bool, optional
            overwrite current state with default values of each parameter, useful for initialization
            or resetting after making changes, by default False
        loud : bool, optional
            flag to print which parameters are being set/updated, by default True

        Keywords
        --------
            use_colorbar: bool, optional
                flag to add a colorbar to the image. not currently implemented for 
                two-color images, by default False
            cbar_label: str, optional
                label for the colorbar, by default ''
            cbar_logspace: bool, optional
                flag for whether ticks on the colorbar should be log-spaced, by default False


        Example usage
        -------------
        ```python 
        gasStudio.set_ImageParams(
            use_colorbar=True,
            cbar_label='Temperature',
            cbar_logspace=True,
            figure_label='t = 13.8 Gyr')```"""

        default_kwargs = {
            'use_colorbar':False,
            'cbar_label':'',
            'cbar_logspace':True,
            }

        for kwarg in list(kwargs):
            ## only set it here if it was passed
            if kwarg in default_kwargs:
                ## remove it from default_kwargs
                default_kwargs.pop(kwarg)
                value = kwargs[kwarg]

                if loud and self.master_loud:
                    print("setting",kwarg,
                        'to user value of:',value)
                ## set it to the object
                setattr(self,kwarg,value)

                ## remove this kwarg to prevent a confusing print message
                ##  suggesting they are being ignored in the parent class
                kwargs.pop(kwarg)
            else:
                if (kwarg not in Studio.set_ImageParams.default_kwargs
                    and kwarg != 'this_setup_id'):
                    if self.master_loud:
                        print(kwarg,'ignored. Did you mean something else?',
                            default_kwargs.keys())

        if use_defaults:
            ## set the remaining image parameters to their default values
            for default_arg in default_kwargs:
                value = default_kwargs[default_arg]
                if loud and self.master_loud:
                    print("setting",default_arg,
                        'to default value of:',value)
                setattr(self,default_arg,value)


        ## set any other image params here
        super().set_ImageParams(use_defaults=use_defaults,loud=loud,**kwargs)

    set_ImageParams.default_kwargs = {
            'use_colorbar':False,
            'cbar_label':'',
            'cbar_logspace':True}
    set_ImageParams.default_kwargs.update(
        Studio.set_ImageParams.default_kwargs)

    append_function_docstring(set_ImageParams,Studio.set_ImageParams,prepend_string='passes `kwargs` to:\n')


    def print_ImageParams(self):
        """ Prints the current image parameters."""

        default_kwargs = {
            'use_colorbar':False,
            'cbar_label':'',
            'cbar_logspace':True,
            }

        ## print the current value, not the default value
        for arg in default_kwargs:
            print(arg,'=',getattr(self,arg))

        ## call the super class' print image params
        super().print_ImageParams()

    def projectAlongLOS(
        self,
        weights:np.ndarray,
        weight_name:str,
        quantities:np.ndarray,
        quantity_name:str,
        use_metadata:bool=True,
        save_meta:bool=True,
        assert_cached:bool=False,
        loud:bool=True,
        **kwargs, 
        ):
        """Projects a weighted quantity along the LOS into pixels. Projection is 
            done with a cubic spline kernel that is renormalized to conserve mass. 

        I.e. 

        `renorm_i = sum_j ( k(r_ij,h_i))`

        where j is a sum over the pixels particle i contributes to.

        The maps computed in the flattened pixel array at index j are then:

        `W[j] = sum_i( k(r_ij,h_i)/renorm_i * weight)`
        `Q[j] = sum_i( k(r_ij,h_i)/renorm_i * weight * quantity) / W[j]`

        Parameters
        ----------
        weights : np.ndarray, optional
            array of weights to use alongside the smoothing kernel
            (typically `self.snapdict['Masses']`).
        weight_name : str
            Name of the key in the snapdict that should be used as the weights
            if weights are not passed.
            special weight_names are `Volumes` and `Ones`, which do not 
            have to be present in the snapdict.
        quantities : np.ndarray, optional
            array of quantities to project along the line of sight.
        quantity_name : str
            Name of the field that is being projected, should be in the snapdict if 
            quantities is not passed.
        use_metadata : bool, optional
            flag for whether a cached result should be used (if it exists), by default True
        save_meta : bool, optional
            flag to save the result in the cache, by default True
        assert_cached : bool, optional
            flag to require a cache hit and raise an exception otherwise, by default False
        loud : bool, optional
            flag for whether cache hits/misses should be announced
            to the console., by default True

        Returns
        -------
        weightMap: np.ndarray
            image of sum(weights) along the line of sight in each pixel
        weightWeightedQuantityMap: np.ndarray
            image of sum(weights*quantity) along the line of sight in each pixel
        """
        

        @metadata_cache(
            self.this_setup_id,  ## hdf5 file group name
            ['%sMap'%weight_name.lower(),
                '%sWeighted%sMap'%(
                weight_name.lower(),
                quantity_name.title())],
            use_metadata=use_metadata,
            save_meta=save_meta,
            assert_cached=assert_cached,
            loud=loud,
            force_from_file=True) ## read from cache file, not attribute of object
        def inner_projectAlongLOS(
            self:GasStudio,
            weights:np.ndarray,
            weight_name:str,
            quantities:np.ndarray,
            quantity_name:str,
            snapdict_name:str='gas'):

            ## make the actual C call
            weightMap, weightWeightedQuantityMap = getImageGrid(
                *self.__prepareCoordinates( ## extract coordinates, weights, quantities from snapdict
                    snapdict_name,
                    weights,
                    weight_name,
                    quantities,
                    quantity_name),
                    loud=self.master_loud)

            if self.master_loud:
                print('-done')

            return weightMap, weightWeightedQuantityMap

        return inner_projectAlongLOS(
            self,
            weights,
            weight_name,
            quantities,
            quantity_name,
            **kwargs)

    def __prepareCoordinates(
        self,
        snapdict_name:str,
        weights:np.ndarray,
        weight_name:str,
        quantities:np.ndarray,
        quantity_name:str):
        """ 

        Parameters
        ----------
        snapdict_name : str
            [description]
        weights : np.ndarray
            array of weights to use alongside the smoothing kernel
            (typically `self.snapdict['Masses']`).
        weight_name : str
            Name of the key in the snapdict that should be used as the weights
            if weights are not passed.
            special weight_names are `Volumes` and `Ones`, which do not 
            have to be present in the snapdict.
        quantities : np.ndarray
            array of quantities to project along the line of sight.
        quantity_name : str
            Name of the field that is being projected, should be in the snapdict if 
            quantities is not passed.

        Returns
        -------
        self.Xmin : float
            minimum x value of coordinates, already applied just returned for convenience
        self.Xmax : float
            maximum x value of coordinates, already applied just returned for convenience
        self.Ymin : float
            minimum y value of coordinates, already applied just returned for convenience
        self.Ymax : float
            maximum y value of coordinates, already applied just returned for convenience
        self.npix_x : float
            number of pixels along the x axis
        self.npix_y : float
            number of pixels along the y axis
        pos : np.ndarray
            coordinate data for particles in kpc 
        weights : np.ndarray
            filtered array of weights to use alongside the smoothing kernel
            (typically `self.snapdict['Masses']`).
        quantities : np.ndarray
            filtered array of quantities to project along the line of sight.
        hsml : np.ndarray
            filtered array of smoothing lengths/radii for each particle

        Raises
        ------
        KeyError
            if `snapdict_name` is not `'gas'` or `'star'`
        """
        

        ## pick which particle type we're projecting, 
        if snapdict_name != 'gas' and snapdict_name != 'star':
            raise KeyError("Choose between 'gas' or 'star' snapdict!")
        
        full_snapdict_name = '%s_snapdict'%snapdict_name
        
        ## use the masked version of the snapdict if it was passed
        if hasattr(self,'masked_'+full_snapdict_name):
            if self.master_loud:
                print("Used masked_snapdict, delete it if you don't want it anymore")
            full_snapdict_name = 'masked_'+full_snapdict_name

        snapdict = getattr(self,full_snapdict_name)

        ## unpack the snapshot data from the snapdict
        Coordinates = snapdict['Coordinates'] ## kpc

        if "SmoothingLength" not in snapdict:
            Hsml = self.get_HSML(snapdict_name)
            if 'masked_' in full_snapdict_name:
                Hsml = Hsml[self.mask]
        else:
            Hsml = snapdict['SmoothingLength'] ## kpc

        ## only important if you are neighbor finding and you want a periodic box.
        ##  for most purposes, you don't. 
        BoxSize = 1000 #snapdict['BoxSize'] 

        if weights is None:
            ## account for possibility of volume weighting
            if weight_name not in snapdict:
                if weight_name == 'Volume':
                    weights = 4/3 * np.pi*Hsml**3 / 32 ## kpc^3
                elif weight_name == 'Ones':
                    weights = np.ones(Coordinates.shape[0])
                else:
                    raise KeyError(weight_name,'is not in gas_snapdict')
            else:
                weights = snapdict[weight_name]

        if quantities is None:
            if quantity_name not in snapdict:
                ## need to rotate the velocities in here
                if quantity_name in ['Vx','Vy','Vz']:
                    vels = self.camera.rotate_array(snapdict['Velocities'])
                    quantities = vels[:,['Vx','Vy','Vz'].index(quantity_name)]
                ## was passed something that we don't know what to do with
                elif 'Log' in quantity_name or 'log' in quantity_name:
                    upper = quantity_name.replace('Log','')
                    lower = quantity_name.replace('log','')
                    if upper in snapdict:
                        quantities = np.log10(snapdict[upper])
                    elif lower in snapdict:
                        quantities = np.log10(snapdict[lower])
                    else:
                        raise KeyError(quantity_name,'is not in gas_snapdict')
                else:
                    raise KeyError(quantity_name,'is not in gas_snapdict')
            else:
                quantities = snapdict[quantity_name]

        ## cull the particles outside the frame and cast to float32
        pos,box_mask = self.camera.project_and_clip(pos)

        if self.master_loud: print("projecting %d particles"%np.sum(box_mask))

        weights = weights[box_mask].astype(np.float32)
        quantities = quantities[box_mask].astype(np.float32)
        hsml = Hsml[box_mask].astype(np.float32)

        return (
            self.Xmin,self.Xmax,
            self.Ymin,self.Ymax,
            self.npix_x,self.npix_y,
            pos,weights,quantities,
            hsml)

    def __quick_projectAlongLOS(
        self,
        weights:np.ndarray,
        weight_name:str,
        quantities:np.ndarray,
        quantity_name:str,
        snapdict_name='gas',
        loud:bool=True,
        **kwargs):
        """ approximates the projection of a weighted quantity along the LOS into pixels.
            using a 2d histogram and assuming particles have 0 extent. much faster than the 
            real thing and useful for quick and dirty testing.

        Parameters
        ----------
        weights : np.ndarray, optional
            array of weights to use alongside the smoothing kernel
            (typically `self.snapdict['Masses']`).
        weight_name : str
            Name of the key in the snapdict that should be used as the weights
            if weights are not passed.
            special weight_names are `Volumes` and `Ones`, which do not 
            have to be present in the snapdict.
        quantities : np.ndarray, optional
            array of quantities to project along the line of sight.
        quantity_name : str
            Name of the field that is being projected, should be in the snapdict if 
            quantities is not passed.
        use_metadata : bool, optional
            flag for whether a cached result should be used (if it exists), by default True
        save_meta : bool, optional
            flag to save the result in the cache, by default True
        assert_cached : bool, optional
            flag to require a cache hit and raise an exception otherwise, by default False
        loud : bool, optional
            flag for whether cache hits/misses should be announced
            to the console., by default True

        Returns
        -------
        weightMap: np.ndarray
            image of sum(weights) along the line of sight in each pixel
        weightWeightedQuantityMap: np.ndarray
            image of sum(weights*quantity) along the line of sight in each pixel
        """

        (Xmin,Xmax,
        Ymin,Ymax,
        npix_x,npix_y,
        prep_pos,prep_weights,prep_quantities,
        prep_hsml) = self.__prepareCoordinates(
            snapdict_name,
            weights,
            weight_name,
            quantities,
            quantity_name)

        xedges = np.linspace(Xmin,Xmax,npix_x,endpoint=True)
        yedges = np.linspace(Ymin,Ymax,npix_y,endpoint=True)

        xs,ys,_ = prep_pos.T
        weightMap,xedges,yedges = np.histogram2d(
            xs,ys,
            bins=[xedges,yedges],
            weights=prep_weights)

        weightWeightedQuantityMap,xedges,yedges = np.histogram2d(
            xs,ys,
            bins=[xedges,yedges],
            weights=prep_weights*prep_quantities)

        return weightMap, weightWeightedQuantityMap/weightMap

####### produceImage implementation #######
    def produceImage(
        self,
        weight_name:str='Masses',
        quantity_name:str='Temperature',
        weights:np.ndarray=None,quantities:np.ndarray=None,
        min_weight:float=None,max_weight:float=None,
        min_quantity:float=None,max_quantity:float=None,
        weight_adjustment_function=None,
        quantity_adjustment_function=None,
        cmap:str='viridis', ## what colormap to use
        quick:bool=False,
        **kwargs
        ):
        """ Generates a projected image using the stored image parameters.
        Specify whether the image should be "one-color" (i.e. a simple projection map)
        or "two-color" (a projection map where the hue is set by the
        LOS projected weighted quantity and the saturation is set by the LOS projected weight)
        by passing values for `min_weight`, `max_weight`, `min_quantity`, and `max_quantity`.

        Parameters
        ----------
        weight_name : str, optional
            Name of the key in the snapdict that should be used as the weights
            if weights are not passed. Special weight_names are `Volumes` and `Ones`,
            which do not have to be present in the snapdict, by default 'Masses'
        quantity_name : str, optional
            Name of the field that is being projected, should be in the snapdict if 
            quantities is not passed, by default 'Temperature'
        weights : np.ndarray, optional
            array of weights to use alongside the smoothing kernel
            (typically `self.snapdict['Masses']`). 
            If None searches `self.snapdict` for `weight_name`, by default None
        quantities : np.ndarray, optional
            array of quantities to project along the line of sight.
            If None searches `self.snapdict` for `weight_name`, by default None
        min_weight : float, optional
            minimum value to use for the colorbar of weight projection maps
            and saturation for two-color images, by default None
        max_weight : float, optional
            maximum value to use for the colorbar of weight projection maps
            and saturation for two-color images, by default None
        min_quantity : float, optional
            minimum value to use for the colorbar of weighted quantity projection maps
            and hue for two-color images, by default None
        max_quantity : float, optional
            maximum value to use for the colorbar of weighted quantity projection maps
            and hue for two-color images, by default None
        weight_adjustment_function : function, optional
            function to apply to the `weightMap` after it is
            returned by `~firestudio.studios.gas_studio.GasStudio.projectAlongLOS`
            for the purposes of specifying `min_weight` and `max_weight`, by default None
        quantity_adjustment_function : function, optional
            function to apply to the `weightWeightedQuantityMap` after it is
            returned by `~firestudio.studios.gas_studio.GasStudio.projectAlongLOS`
            for the purposes of specifying `min_quantity` and `max_quantity`, by default None
        cmap : str, optional
            name of colormap to apply to image, by default 'viridis'

        Returns
        -------
        np.ndarray(npix_x,npix_y,3)
            array of RGB image pixel values

        Raises
        ------
        ValueError
            if an invalid combination of min/max weight/quantity parameters are 
            passed.

        Example usage
        -------------
        ```python
        ## makes a gas surface density map
        gasStudio.render(
            weight_name='Masses',
            min_weight=-0.1,
            max_weight=1.5,
            weight_adjustment_function= lambda x: np.log10(x/gasStudio.Acell)+10-6 ## log10(msun/pc^2)
            )

        ## makes a mass weighted temperature map
        gasStudio.render(
            weight_name='Masses',
            quantity_name='Temperature',
            min_quantity=2,
            max_quantity=7,
            quantity_adjustment_function= np.log10
            )

        ## makes a saturation-hue gas surface density + Temperature map
        gasStudio.render(
            weight_name='Masses',
            min_weight=-0.1,
            max_weight=1.5,
            weight_adjustment_function= lambda x: np.log10(x/gasStudio.Acell)+10-6 ## log10(msun/pc^2)
            quantity_name='Temperature',
            min_quantity=2,
            max_quantity=7,
            quantity_adjustment_function= np.log10
            )```"""

        self.cmap = cmap


        ## load the requested maps
        if not quick:
            weightMap, weightWeightedQuantityMap = self.projectAlongLOS(
                weights,
                weight_name,
                quantities,
                quantity_name,
                **kwargs)
        else:
            weightMap, weightWeightedQuantityMap = self.__quick_projectAlongLOS(
                weights,
                weight_name,
                quantities,
                quantity_name,
                **kwargs)

        ## apply any unit corrections, take logs, etc...
        if weight_adjustment_function is not None:
            weightMap = weight_adjustment_function(weightMap)

        if quantity_adjustment_function is not None: 
            weightWeightedQuantityMap = quantity_adjustment_function(weightWeightedQuantityMap)

        ## plot a hue-brightness image, convert to 0->1 space
        if (min_weight is not None and 
            max_weight is not None and 
            min_quantity is not None and
            max_quantity is not None):

            image_W = self.renormalizeTransposeImage(
                weightMap, 
                min_weight,max_weight,
                weight_name)

            image_Q = self.renormalizeTransposeImage(
                weightWeightedQuantityMap,
                min_quantity,max_quantity,
                quantity_name)

            #self.cbar_label = 'ERROR'
            self.cbar_min = min_quantity
            self.cbar_max = max_quantity

            if self.master_loud:
                print("TODO:Need to create a 2-axis colorbar.")

        ## plot a weight map, convert to 0->1 space
        elif (min_weight is not None and 
            max_weight is not None):

            #self.cbar_label = 'los %s ' % (
                #self.weight_name,
                #self.quantity_name.title())

            image_Q = self.renormalizeTransposeImage(
                weightMap, 
                min_weight,max_weight,
                weight_name)

            image_W = None

            self.cbar_min = min_weight
            self.cbar_max = max_weight
        
        ## plot a quantity map, convert to 0->1 space
        elif (min_quantity is not None and
            max_quantity is not None):

            #self.cbar_label = 'los %s-weighted %s' % (
                #self.weight_name,
                #self.quantity_name.title())

            image_Q = self.renormalizeTransposeImage(
                weightWeightedQuantityMap,
                min_quantity,max_quantity,
                quantity_name)

            image_W = None

            self.cbar_min = min_quantity
            self.cbar_max = max_quantity

        else: raise ValueError("Use (min/max)_(weight/quantity) kwargs to set image")

        ## convert the images from 0->1 space to 0-> 255 space
        final_image = produce_cmap_hsv_image(image_Q, image_W, cmap=self.cmap) 

        return final_image

####### plotImage implementation #######
    def plotImage(
        self,
        ax:plt.Axes,
        final_image:np.ndarray):
        """Implementation of overlaying artists on top of projected image.
        if `self.use_colorbar==True`
            we add a colorbar to the image
        See `~firestudio.studios.gas_studio.GasStudio.set_ImageParams` for details.

        Also calls the base method `firestudio.studios.studio.Studio.plotImage`. 

        Parameters
        ----------
        ax : plt.Axes
            axis to plot image to 
        final_image : np.ndarray
            array of RGB image pixel values
        """

        ## run Studio's plotImage method
        super().plotImage(ax,final_image)

        ## colour bar
        if self.use_colorbar:
            ## do we need to exponentiate the cbar limits?
            if self.cbar_logspace:
                cbar_min,cbar_max = 10**self.cbar_min,10**self.cbar_max 
            else:
                cbar_min,cbar_max = self.cbar_min,self.cbar_max 
        
            addColorbar(
                ax,get_cmap(self.cmap),
                cbar_min,cbar_max,
                self.cbar_label,
                logflag = self.cbar_logspace,
                fontsize=self.fontsize,
                cmap_number=0.25)

def getImageGrid(
    Xmin,Xmax,
    Ymin,Ymax,
    npix_x,npix_y,
    pos,weight,quantity,
    hsml,
    loud=True):

    ## set c-routine variables
    n_smooth = pos.shape[0]

    ## output array for sum along the line of sight
    weightMap = np.zeros(shape = (npix_x,npix_y),dtype=np.float32)

    ## output array for average along the line of sight
    weightWeightedQuantityMap = np.zeros(shape = (npix_x,npix_y),dtype=np.float32)
    
    ## create hsml output array
    if hsml is None: raise ValueError(
            "HSML cannot be None and weights != masses." + 
            " We don't check if weights == masses, so we'll just assume" +
            " they're not for ultimate safety.")
    #else: print("Using provided smoothing lengths")
    
    ## make sure everything is in single precision lest we
    ##  make a swiss-cheese magenta nightmare, #neverforget 6/15/17
    c_f_p      = ctypes.POINTER(ctypes.c_float)

    if loud: print('------------------------------------------')
    curpath = os.path.realpath(__file__)
    curpath = os.path.split(curpath)[0] #split off this filename
    curpath = os.path.split(curpath)[0] #split off studios direcotry
    c_obj_path = os.path.join(
        curpath,
        'utils',
        'C_routines',
        'HsmlAndProject_cubicSpline/hsml_project.so')

    if not os.path.isfile(c_obj_path): raise IOError(
        'Missing ' + c_obj_path + ' - compile the missing file and restart.')

    c_obj = ctypes.CDLL(c_obj_path)

    c_obj.hsml_project( 
        ctypes.c_int(n_smooth), ## number of particles
        vfloat(copy.copy(pos[:,0])),vfloat(copy.copy(pos[:,1])), ## x-y positions of particles
        vfloat(hsml),  ## smoothing lengths of star + gas particles
        vfloat(weight), ## attenuation masses of star + gas particles, stars are 0 
        ## emission in each band of star+gas particles, gas is 0 
        vfloat(quantity),  
        ## x-y limits of the image
        ctypes.c_float(Xmin),ctypes.c_float(Xmax),
        ctypes.c_float(Ymin),ctypes.c_float(Ymax), 
        ctypes.c_int(npix_x),ctypes.c_int(npix_y), ## output shape
        weightMap.ctypes.data_as(c_f_p), ## mass map
        weightWeightedQuantityMap.ctypes.data_as(c_f_p))
    
    # convert into Msun/pc^2
    #unitmass_in_g = 1.9890000e+43 
    #solar_mass    = 1.9890000e+33
    #conv_fac = (unitmass_in_g/solar_mass) / (1.0e3)**2 ## Msun/pc^2
    #columnDensityMap *= conv_fac
    if loud:
        print('------------------------------------------')

        print('minmax(weightMap)',
            np.min(weightMap),
            np.max(weightMap))
        print('Fraction deposited:',np.sum(weightMap)/np.sum(weight))

        print('minmax(weightWeightedQuantityMap)',
            np.min(weightWeightedQuantityMap),
            np.min(weightWeightedQuantityMap))
   
    weightWeightedQuantityMap = weightWeightedQuantityMap/weightMap

    return weightMap,weightWeightedQuantityMap

###### Color helper functions
from matplotlib.colors import rgb_to_hsv, hsv_to_rgb 
def produce_colmap(cmap_name):
    cmap = get_cmap(cmap_name)
    ## discretize the colormap into 256 parts...
    return [list(cmap(i/255.)[:3]) for i in range(0,256)]

def produce_cmap_hsv_image(image_1, image_2,cmap='viridis'): 
    # image_1 and image_2 are arrays of pixels 
    # with integer values in the range 0 to 255. 
    # These will be mapped onto hue and brightness, 
    # respectively, with saturation fixed at 1. 
    cols = produce_colmap(cmap)
    cols = np.array(cols) 

    cols_2d = np.ones((len(cols), 1, 3)) 
    cols_2d[:, 0, 0] = cols[:, 0] 
    cols_2d[:, 0, 1] = cols[:, 1] 
    cols_2d[:, 0, 2] = cols[:, 2] 
    cols_hsv = rgb_to_hsv(cols_2d) 
    hue_viridis = cols_hsv[:, 0, 0] 

    npix_x = len(image_1) 
    npix_y = len(image_1[0]) 
    if image_2 is not None:
        output_image_hsv = np.zeros((npix_x, npix_y, 3)) 
        for i in range(npix_x): 
            for j in range(npix_y): 
                output_image_hsv[i, j, 0] = hue_viridis[image_1[i, j]] 
                output_image_hsv[i, j, 1] = 1.0 
                output_image_hsv[i, j, 2] = float(image_2[i, j]) / 255.0 

        output_image_rgb = hsv_to_rgb(output_image_hsv) 
    else:
        output_image_rgb = np.zeros((npix_x, npix_y, 3)) 
        for i in range(npix_x): 
            for j in range(npix_y): 
                output_image_rgb[i,j]=cols[image_1[i,j]]
                
    return output_image_rgb 

