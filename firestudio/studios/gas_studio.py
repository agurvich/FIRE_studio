## builtin imports
import os
import sys 
import h5py
import matplotlib 
matplotlib.use('Agg') 
import numpy as np 
import ctypes

## abg_python imports
from abg_python.all_utils import filterDictionary,append_function_docstring
from abg_python.plot_utils import addColorbar
from abg_python.galaxy.metadata_utils import metadata_cache

## firestudio imports
import firestudio.utils.gas_utils.my_colour_maps as mcm 
from firestudio.studios.studio import Studio

class GasStudio(Studio):
    """ FIREstudio class for making gas projection images.
        Can be used for stars, but you will either have to pass smoothing lengths
        in or allow FIREstudio to calculate them itself, which  can take a long time. 

        Methods:
    """

    def __repr__(self):
        return 'GasStudio instance'

####### makeOutputDirectories implementation #######
    def makeOutputDirectories(self,datadir):
        raise NotImplementedError("rethinking I/O")

        print('Drawing %s:%d'%(self.snapdir,self.snapnum)+' to:%s'%self.datadir)

        ## make the path to store all plots in the datadir
        path = os.path.join(datadir,'Plots')
        if not os.path.isdir(path):
            os.mkdir(path)

        ## make the path to store final images
        self.image_dir = os.path.join(path,'GasTwoColour')
        if not os.path.isdir(self.image_dir):
            os.mkdir(self.image_dir)

        ## make the paths to store projection hdf5 files
        self.projection_dir = os.path.join(path,'Projections')
        if not os.path.isdir(self.projection_dir):
            os.mkdir(self.projection_dir)

    def set_ImageParams(
        self,
        use_defaults=False,
        loud=True,
        **kwargs):
        """ 
            Input: 
                use_colorbar=False
                cbar_label=''
                cbar_logspace=False,
                cbar_max=7"""

        default_kwargs = {
            'use_colorbar':False,
            'cbar_label':'',
            'cbar_logspace':True,
            }

        for kwarg in kwargs:
            ## only set it here if it was passed
            if kwarg in default_kwargs:
                ## remove it from default_kwargs
                default_kwargs.pop(kwarg)
                value = kwargs[kwarg]

                if loud:
                    print("setting",kwarg,
                        'to user value of:',value)
                ## set it to the object
                setattr(self,kwarg,)

        if use_defaults:
            ## set the remaining image parameters to their default values
            for default_arg in default_kwargs:
                value = default_kwargs[default_arg]
                if loud:
                    print("setting",default_arg,
                        'to default value of:',value)
                setattr(self,default_arg,value)

        ## set any other image params here
        super().set_ImageParams(use_defaults=use_defaults,loud=loud,**kwargs)

    append_function_docstring(set_ImageParams,Studio.set_ImageParams)


    def print_ImageParams():
        """ Prints current image setup to console.

            Input:
                None

            Output:

                None """

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


    def weightAvgAlongLOS(
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
        """ 
            Input:
                weights
                weight_name
                quantities
                quantity_name
                use_metadata=True
                save_meta=True
                assert_cached=False
                loud=True
                
                snapdict_name - 'gas' or 'star', where to load arrays that hasn't 
                    been passed from

            Output:
                weightMap - 
                weightWeightedQuantityMap -"""

        @metadata_cache(
            self.this_setup_id,  ## hdf5 file group name
            ['%sMap'%weight_name.lower(),
                '%sWeighted%sMap'%(
                weight_name.lower(),
                quantity_name.title())][:1+(quantity_name!='Ones')],
            use_metadata=use_metadata,
            save_meta=save_meta,
            assert_cached=assert_cached,
            loud=loud)
        def inner_weight_along_los(
            self,
            weights,
            weight_name,
            quantities,
            quantity_name,
            snapdict_name='gas'):

            ## pick which particle type we're projecting, 
            if snapdict_name != 'gas' and snapdict_name != 'star':
                raise ValueError("Choose between 'gas' or 'star' snapdict!")
            
            full_snapdict_name = '%s_snapdict'%snapdict_name
            
            ## use the masked version of the snapdict if it was passed
            if hasattr(self,'masked_'+full_snapdict_name):
                print("Used masked_snapdict, delete it if you don't want it anymore")
                full_snapdict_name = 'masked_'+full_snapdict_name

            snapdict = getattr(self,full_snapdict_name)

            ## unpack the snapshot data from the snapdict
            Coordinates = snapdict['Coordinates'] ## kpc
            
            ## TODO use phil's get star smoothing length thing to 
            ##  compute smoothing lengths of particles (and cache them) 
            ##  if they're missing.
            if "SmoothingLength" not in snapdict:
                Hsml = self.get_HSML(snapdict_name)
                assert type(Hsml) == np.ndarray
                if 'masked_' in full_snapdict_name:
                    Hsml = Hsml[self.mask]
            else:
                Hsml = snapdict['SmoothingLength'] ## kpc

            ## only important if you are neighbor finding and you want a periodic box.
            ##  for most purposes, you don't. 
            BoxSize = snapdict['BoxSize'] 

            if weights is None:
                ## account for possibility of volume weighting
                if weight_name not in snapdict:
                    if weight_name == 'Volume':
                        weights = 4/3 * np.pi*Hsml**3 / 32 ## kpc^3
                    else:
                        raise KeyError(weight_name,'is not in gas_snapdict')
                else:
                    weights = snapdict[weight_name]

            if quantities is None:
                if quantity_name not in snapdict:
                    if quantity_name == 'Ones':
                        quantities = np.ones(weights.size)
                    else:
                        raise KeyError(quantity_name,'is not in gas_snapdict')

                else:
                    quantities = snapdict[quantity_name]

            
            ## cull the particles outside the frame and cast to float32
            box_mask = self.cullFrameIndices(Coordinates) ## TODO is this where I want to rotate?

            print("projecting %d particles"%np.sum(box_mask))

            pos = Coordinates[box_mask].astype(np.float32)
            weights = weights[box_mask].astype(np.float32)
            quantities = quantities[box_mask].astype(np.float32)
            hsml = Hsml[box_mask].astype(np.float32)

            frame_center = self.frame_center.astype(np.float32)

            ## rotate by euler angles if necessary
            pos = self.rotateEuler(self.theta,self.phi,self.psi,pos)



            ## make the actual C call
            weightMap, weightWeightedQuantityMap = getImageGrid(
                BoxSize,
                self.Xmin,self.Xmax,
                self.Ymin,self.Ymax,
                self.Zmin,self.Zmax,
                self.npix_x,self.npix_y,
                pos,weights,quantities,
                hsml)

            print('-done')

            return_list = [weightMap, weightWeightedQuantityMap] ## lol

            return return_list[:1+(quantity_name!='Ones')]

        return_value = inner_weight_along_los(
            self,
            weights,
            weight_name,
            quantities,
            quantity_name,
            **kwargs)

        if quantity_name == 'Ones':
            return return_value[0],return_value[0] ## return the same map twice
        else:
            return return_value

    def volumeWeightAlongLOS(
        self,
        quantity,
        quantity_name, 
        **kwargs):
        """ Wrapper function for easier API if quantity = None 
            but quantity_name is in the snapdict then it will be read for you.

            Volume is calculated as 4/3 pi hsml^3 / 32.
            If you'd like a different volume weight then set self.gas_snapdict['Volume'] = volumes.

            Alternatively, call weightAvgAlongLOS directly with quantity and weights passed.

            Input:
                quantity
                quantity_name

                ---
                
                kwargs passed:
                    use_metadata=True
                    save_meta=True
                    assert_cached=False
                    loud=True
                    snapdict_name='gas'

            Output: 
                volumeMap - a map of the deposited cell volume along the LOS
                    in each pixel.
                volumeWeightedQuantityMap - the volume weighted quantity along the
                    LOS in each pixel."""

        return self.weightAvgAlongLOS(
            None, ## read it from the snapdict
            'Volume',
            quantity,
            quantity_name,
            **kwargs)

    def massWeightAlongLOS(
        self,
        quantity,
        quantity_name,
        **kwargs):
        """ Wrapper function for easier API if quantity = None 
            but quantity_name is in the snapdict then it will be read for you.

            Input:
                quantity
                quantity_name

                kwargs passed:
                    use_metadata=True
                    save_meta=True
                    assert_cached=False
                    loud=True
                    snapdict_name='gas'
            
            Output:
                massMap"""

        return self.weightAvgAlongLOS(
            None, ## read it from the snapdict.
            'Masses',
            quantity,
            quantity_name,
            **kwargs)

####### produceImage implementation #######
    def produceImage(
        self,
        weight_name='Masses',
        quantity_name='Temperature',
        weights=None,quantities=None,
        min_weight=None,max_weight=None,
        min_quantity=None,max_quantity=None,
        weight_adjustment_function=None,
        quantity_adjustment_function=None,
        use_colorbar=False,
        cmap='viridis', ## what colormap to use
        **kwargs
        ):

        self.cmap = cmap

        ## load the requested maps
        weightMap, weightWeightedQuantityMap = self.weightAvgAlongLOS(
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

        else:
            raise ValueError("Use (min/max)_(weight/quantity) kwargs to set image")


        ## convert the images from 0->1 space to 0-> 255 space
        final_image = mcm.produce_cmap_hsv_image(image_Q, image_W, cmap=self.cmap) 

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
                ax,mcm.get_cmap(self.cmap),
                cbar_min,cbar_max,
                self.cbar_label,
                logflag = self.cbar_logspace,
                fontsize=self.fontsize,
                cmap_number=0.25)

append_function_docstring(GasStudio,GasStudio.massWeightAlongLOS)
append_function_docstring(GasStudio,GasStudio.volumeWeightAlongLOS)
append_function_docstring(GasStudio,GasStudio.weightAvgAlongLOS)
append_function_docstring(GasStudio,GasStudio.set_ImageParams)
append_function_docstring(GasStudio,Studio)

def getImageGrid(
    BoxSize,
    Xmin,Xmax,
    Ymin,Ymax,
    Zmin,Zmax,
    npix_x,npix_y,
    pos,mass,quantity,
    hsml):

    ## set c-routine variables
    desngb   = 32
    Axis1    = 0
    Axis2    = 1
    Axis3    = 2

    Hmax     = 0.5*(Xmax-Xmin) ## ignored if smoothing lengths are passed in

    n_smooth = pos.shape[0]

    ## output array for sum along the line of sight
    weightMap = np.zeros(shape = (npix_x,npix_y),dtype=np.float32)

    ## output array for average along the line of sight
    weightWeightedQuantityMap = np.zeros(shape = (npix_x,npix_y),dtype=np.float32)
    
    ## create hsml output array
    if hsml is None:
        raise ValueError("HSML cannot be None and weights != masses.",
            "We don't check if weights == masses, so we'll just assume",
            "they're not for ultimate safety.")
        hsml = np.zeros(mass.shape[0],dtype=np.float32)
    #else:
        #print("Using provided smoothing lengths")
    
    ## make sure everything is in single precision lest we
    ##  make a swiss-cheese magenta nightmare, #neverforget 6/15/17
    c_f_p      = ctypes.POINTER(ctypes.c_float)
    pos_p      = pos.astype(np.float32).ctypes.data_as(c_f_p)
    hsml_p     = hsml.astype(np.float32).ctypes.data_as(c_f_p)
    mass_p     = mass.astype(np.float32).ctypes.data_as(c_f_p)
    quantity_p = quantity.astype(np.float32).ctypes.data_as(c_f_p)

    w_f_p    = weightMap.ctypes.data_as(c_f_p)
    q_f_p    = weightWeightedQuantityMap.ctypes.data_as(c_f_p)

    print('------------------------------------------')
    curpath = os.path.realpath(__file__)
    curpath = os.path.split(curpath)[0] #split off this filename
    curpath = os.path.split(curpath)[0] #split off studios direcotry
    c_obj = ctypes.CDLL(os.path.join(
        curpath,'utils','gas_utils','HsmlAndProject_cubicSpline/HsmlAndProject.so'))

    #print(n_smooth)
    #print(pos_p)
    #print(hsml_p)
    #print(mass_p)
    #print(quantity_p)
    #print(Xmin,Xmax)
    #print(Ymin,Ymax)
    #print(Zmin,Zmax)
    #print(npix_x,npix_y)
    #print(desngb)
    #print(Axis1,Axis2,Axis3)
    #print(Hmax,BoxSize)

    c_obj.findHsmlAndProject(
	ctypes.c_int(n_smooth), ## number of particles
	pos_p,hsml_p,mass_p,quantity_p, ## position, mass, and "quantity" of particles
        ctypes.c_float(Xmin.astype(np.float32)),ctypes.c_float(Xmax.astype(np.float32)), ## xmin/xmax
	ctypes.c_float(Ymin.astype(np.float32)),ctypes.c_float(Ymax.astype(np.float32)), ## ymin/ymax
	ctypes.c_float(Zmin.astype(np.float32)),ctypes.c_float(Zmax.astype(np.float32)), ## zmin/zmax
        ctypes.c_int(npix_x),ctypes.c_int(npix_y), ## npixels
	ctypes.c_int(desngb), ## neighbor depth
        ctypes.c_int(Axis1),ctypes.c_int(Axis2),ctypes.c_int(Axis3), ## axes...?
	ctypes.c_float(Hmax),ctypes.c_double(BoxSize), ## maximum smoothing length and size of box
	w_f_p,q_f_p) ## pointers to output cell-mass and cell-mass-weighted-quantity
    print('------------------------------------------')
    
    # convert into Msun/pc^2
    #unitmass_in_g = 1.9890000e+43 
    #solar_mass    = 1.9890000e+33
    #conv_fac = (unitmass_in_g/solar_mass) / (1.0e3)**2 ## Msun/pc^2
    #columnDensityMap *= conv_fac
    print('minmax(weightMap)',
	np.min(weightMap),
	np.max(weightMap))

    # weightWeightedQuantityMap contains the mass-weighted quantity
    print('minmax(weightWeightedQuantityMap)',
	np.min(weightWeightedQuantityMap),
	np.min(weightWeightedQuantityMap))
   
    return weightMap,weightWeightedQuantityMap
