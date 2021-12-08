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
        use_defaults=False,
        loud=True,
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
                , by default False
            cbar_label: str, optional
                , by default ''
            cbar_logspace: bool, optional
                , by default False


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
        weights,
        weight_name,
        quantities,
        quantity_name,
        use_metadata=True,
        save_meta=True,
        assert_cached=False,
        loud=True,
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
            [description]
        weight_name : str
            Name of the key in the snapdict that should be used as the weights
            if weights are not passed.
            special weight_names are `Volumes` and `Ones`, which do not 
            have to be present in the snapdict.
        quantities : np.ndarray, optional
            [description]
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
            [description]
        weightWeightedQuantityMap: np.ndarray
            [description]
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
            self,
            weights,
            weight_name,
            quantities,
            quantity_name,
            snapdict_name='gas'):

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
        snapdict_name,
        weights,
        weight_name,
        quantities,
        quantity_name):
        """ keys required to function:
            Coordinates
            SmoothingLength (optional, will compute though otherwise)
            weight_name (unless it's Ones or passed in as weights)
            quantity_name (unless it's passed in as quantities)

        Parameters
        ----------
        snapdict_name : [type]
            [description]
        weights : [type]
            [description]
        weight_name : [type]
            [description]
        quantities : [type]
            [description]
        quantity_name : [type]
            [description]

        Returns
        -------
        [type]
            [description]

        Raises
        ------
        KeyError
            [description]
        KeyError
            [description]
        KeyError
            [description]
        KeyError
            [description]
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
            assert type(Hsml) == np.ndarray
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

        ## rotate by euler angles if necessary
        pos = self.camera.rotate_array(Coordinates,offset=True)

        ## cull the particles outside the frame and cast to float32
        box_mask = self.cullFrameIndices(pos)

        if self.master_loud:
            print("projecting %d particles"%np.sum(box_mask))

        pos = pos[box_mask].astype(np.float32)
        weights = weights[box_mask].astype(np.float32)
        quantities = quantities[box_mask].astype(np.float32)
        hsml = Hsml[box_mask].astype(np.float32)

        return (
            BoxSize,
            self.Xmin,self.Xmax,
            self.Ymin,self.Ymax,
            self.Zmin,self.Zmax,
            self.npix_x,self.npix_y,
            pos,weights,quantities,
            hsml)

    def __quick_projectAlongLOS(
        self,
        weights,
        weight_name,
        quantities,
        quantity_name,
        snapdict_name='gas',
        **kwargs):
        """ approximates the projection of a weighted quantity along the LOS into pixels.
            using a 2d histogram and assuming particles have 0 extent. much faster than the 
            real thing and useful for quick and dirty testing.

        Parameters
        ----------
        weights : np.ndarray, optional
            [description]
        weight_name : str
            Name of the key in the snapdict that should be used as the weights
            if weights are not passed.
            special weight_names are `Volumes` and `Ones`, which do not 
            have to be present in the snapdict.
        quantities : np.ndarray, optional
            [description]
        quantity_name : str
            Name of the field that is being projected, should be in the snapdict if 
            quantities is not passed.
        snapdict_name : str, optional
            [description], by default 'gas'

        Returns
        -------
        weightMap: np.ndarray
            [description]
        weightWeightedQuantityMap: np.ndarray
            [description]
        """

        (BoxSize,
        Xmin,Xmax,
        Ymin,Ymax,
        Zmin,Zmax,
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
    def render(
        self,
        ax=None,
        **kwargs):
        """ Plots a projected image using the stored image parameters.

        Parameters
        ----------
        ax : matplotlib axis, optional
            [description], by default None

        Keywords
        --------
            ax = None -- axis to plot image to, if None will create a new figure

            weight_name = 'Masses' --
            quantity_name = 'Temperature' --
            weights = None -- 
            quantities = None --
            min_weight = None,max_weight = None --
            min_quantity = None,max_quantity = None --
            weight_adjustment_function = None --
            quantity_adjustment_function = None --
            use_colorbar = False --
            cmap = 'viridis' -- what colormap to use
            quick = False -- flag to use a simple 2d histogram (for comparison or
                for quick iteration as the user defines the image parameters)

        Returns
        -------
            ax -- the axis the image was plotted to
            final_image -- 2x2x3 RGB pixel array

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

        if ax is None: fig,ax = plt.figure(),plt.gca()
        else: fig = ax.get_figure()

        ## remap the C output to RGB space
        final_image = self.produceImage(**kwargs)

        ## plot that RGB image and overlay scale bars/text
        self.plotImage(ax,final_image)

        ## save the image
        if self.savefig is not None:
            self.saveFigure(ax,self.savefig)

        return ax,final_image

    def produceImage(
        self,
        weight_name='Masses',
        quantity_name='Temperature',
        weights=None,quantities=None,
        min_weight=None,max_weight=None,
        min_quantity=None,max_quantity=None,
        weight_adjustment_function=None,
        quantity_adjustment_function=None,
        cmap='viridis', ## what colormap to use
        quick=False,
        **kwargs
        ):
        """[summary]

        Parameters
        ----------
        weight_name : str, optional
            [description], by default 'Masses'
        quantity_name : str, optional
            [description], by default 'Temperature'
        weights : [type], optional
            [description], by default None
        quantities : [type], optional
            [description], by default None
        min_weight : [type], optional
            [description], by default None
        max_weight : [type], optional
            [description], by default None
        min_quantity : [type], optional
            [description], by default None
        max_quantity : [type], optional
            [description], by default None
        weight_adjustment_function : [type], optional
            [description], by default None
        quantity_adjustment_function : [type], optional
            [description], by default None
        use_colorbar : bool, optional
            [description], by default False
        cmap : str, optional
            [description], by default 'viridis'

        Returns
        -------
        [type]
            [description]

        Raises
        ------
        ValueError
            [description]
        """

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
        ax,
        final_image
        ):

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
    BoxSize,
    Xmin,Xmax,
    Ymin,Ymax,
    Zmin,Zmax,
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

