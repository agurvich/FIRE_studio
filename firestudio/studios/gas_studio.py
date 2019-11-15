## builtin imports
import os
import sys 
import h5py
import matplotlib 
matplotlib.use('Agg') 
import numpy as np 
import ctypes

## abg_python imports
from abg_python.all_utils import filterDictionary
from abg_python.plot_utils import addColorbar

## firestudio imports
import firestudio.utils.gas_utils.my_colour_maps as mcm 
from firestudio.studios.studio import Studio

class GasStudio(Studio):
    """
    Input:
        snapdir - location that the snapshots live in
        snapnum - snapshot number
        datadir - directory to output the the intermediate grid files and output png
        frame_half_width - half-width of image in data space
        frame_depth - half-depth of image in data space 

    Optional:
        min_den=-0.4 - the minimum of the density color scale
        max_den=1.6 - the maximum of the density color scale
        min_quantity=2 - the minimum of the temperature color scale
        max_quantity=7 - the maximum of the temperature color scale

        single_image=None - string, if it's "Density" it will plot a column 
            density projection, if it's anything else it will be a mass weighted
            `quantity_name` projection. None will be a "two-colour" projection
            with hue determined by `quantity_name` and saturation by density

        quantity_name='Temperature' - the name of the quantity that you're mass weighting
            should match whatever array you're passing in as quantity

        cmap='viridis' - string name for cmap to use 
        use_colorbar=False - flag for whether to plot a colorbar at all
        cbar_label=None - string that shows up next to colorbar, good for specifying 
            units, otherwise will default to just quantity_name.title()
        take_log_of_quantity=True - should we plot the log of the resulting quantity map?
        Hsml=None - Provided, gas smoothing lengths, speeds up projection calculation
        snapdict=None - Dictionary-like holding gas snapshot data, open from disk if None
        use_hsml=True - Flag to use the provided Hsml argument (implemented to test speedup)
        intermediate_file_name = "proj_maps" ##  the name of the file to save maps to
    """ + "------- Studio\n" + Studio.__doc__

    def __init__(
        self,
        snapdir,snapnum, ## snapshot directory and snapshot number
        datadir, ## directory to put intermediate and output files
        frame_half_width, ## half-width of image in x direction
        frame_depth, ## z-depth of image (thickness is 2*frame_depth)
        min_den=-0.4, ## the minimum of the density color/saturation scale
        max_den=1.6, ## the maximum of the density color/saturation scale
        min_quantity=2, ## the minimum of the temperature color scale
        max_quantity=7, ## the maximum of the temperature color scale
        single_image = None, ## string to determine what sort of 1-color image to make
        quantity_name='Temperature', ## quantity to make a mass weighted map/2 color image
        take_log_of_quantity=True, ## take log of mass weighted quantity map?
        cmap='viridis', ## what colormap to use
        use_colorbar = False,
        use_hsml = True, ## flag to use the smoothing lengths passed
        snapdict = None, ## provide an open snapshot dictionary to save time opening
        **kwargs):

        ## image limits
        self.min_den,self.max_den = min_den,max_den
        self.min_quantity,self.max_quantity = min_quantity,max_quantity

        ## what is the quantity we want to make a 2 color image with?
        ##  (or 1 color mass weighted map)
        self.quantity_name=quantity_name
        self.take_log_of_quantity=take_log_of_quantity
        self.single_image = single_image
        self.cmap = cmap
        self.use_colorbar = use_colorbar

        self.use_hsml = use_hsml

        ## call Studio's init
        super().__init__(
            snapdir,snapnum,
            datadir,
            frame_half_width,
            frame_depth,
            **kwargs)

        self.snapdict = snapdict

####### makeOutputDirectories implementation #######
    def makeOutputDirectories(self,datadir):
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

####### projectImage implementation #######
    def projectImage(self,image_names):

        ## open snapshot data if necessary
        if self.snapdict is None:
            self.openSnapshot(
                keys_to_extract = 
                    ['Coordinates',
                    'Masses',
                    'Velocities',
                    self.quantity_name,
                    'SmoothingLength'])

        ## unpack the snapshot data from the snapdict
        Coordinates = self.snapdict['Coordinates']
        Masses = self.snapdict['Masses']
        Quantity = self.snapdict[self.quantity_name]

        if 'SmoothingLength' in self.snapdict and self.use_hsml:
            Hsml = self.snapdict['SmoothingLength']
        else:
            Hsml = None

        BoxSize = self.snapdict['BoxSize']

        ## cull the particles outside the frame and cast to float32
        ind_box = self.cullFrameIndices(Coordinates)

        pos = Coordinates[ind_box].astype(np.float32)
        mass = Masses[ind_box].astype(np.float32)
        quantity = Quantity[ind_box].astype(np.float32)
        hsml = Hsml[ind_box].astype(np.float32) if Hsml is not None else Hsml

        frame_center = self.frame_center.astype(np.float32)

        print('-done')

        ## rotate by euler angles if necessary
        pos = self.rotateEuler(self.theta,self.phi,self.psi,pos)

        ## make the actual C call
        columnDensityMap, massWeightedQuantityMap = getImageGrid(
            BoxSize,
            self.Xmin,self.Xmax,
            self.Ymin,self.Ymax,
            self.Zmin,self.Zmax,
            self.npix_x,self.npix_y,
            pos,mass,quantity,
            self.take_log_of_quantity,
            hsml = hsml)

        ## write the output to an .hdf5 file
        self.writeImageGrid(
            columnDensityMap,
            'columnDensityMap',
            overwrite=self.overwrite)

        self.writeImageGrid(
            massWeightedQuantityMap, 
            "massWeighted%sMap"%self.quantity_name.title(),
            overwrite=self.overwrite)

####### produceImage implementation #######
    def produceImage(self,image_names):
        ## open the hdf5 file and load the maps
        with h5py.File(self.projection_file, "r") as handle:
            this_group=handle[self.this_setup_id]
            columnDensityMap = np.array(this_group['columnDensityMap'])
            massWeightedQuantityMap = np.array(this_group['massWeighted%sMap'%self.quantity_name.title()])

        ## make sure that the maps we're loading are the correct shape
        try:
            assert columnDensityMap.shape == (self.npix_x,self.npix_y)
            assert massWeightedQuantityMap.shape == (self.npix_x,self.npix_y)
        except:
            raise ValueError("Map (%d,%d) is not the correct shape (%d,%d)"%
                (columnDensityMap.shape[0],columnDensityMap.shape[1],self.npix_x,self.npix_y))


        image_rho = self.renormalizeTransposeImage(
            columnDensityMap,
            self.min_den,self.max_den,
            'rho')

        ## only load and renormalize the image_Q if necessary
        if self.single_image != 'Density':
            image_Q = self.renormalizeTransposeImage(
                massWeightedQuantityMap,
                self.min_quantity,self.max_quantity,
                self.quantity_name)

        if self.single_image is None:
            ## Now take the rho and T images, and combine them 
            ##	to produce the final image array. 
            final_image = mcm.produce_cmap_hsv_image(image_Q, image_rho,cmap=self.cmap) 
            self.cbar_label = 'ERROR'

        ## make a column density map
        elif self.single_image == 'Density':
            final_image = mcm.produce_cmap_hsv_image(image_rho,None,cmap=self.cmap)
            self.cbar_label='Column Density (M$_\odot$/pc$^2$)'
            ## set the quantity limits to be the density limits for the colorbar...
            self.min_quantity = self.min_den
            self.max_quantity = self.max_den
            self.take_log_of_quantity=True
        ## make a mass weighted quantity map
        else:
            final_image = mcm.produce_cmap_hsv_image(image_Q,None,cmap=self.cmap)
            self.cbar_label = self.quantity_name.title()

        self.final_image = final_image

        return final_image

####### plotImage implementation #######
    def plotImage(self,ax,image_names):
        ## run Studio's plotImage method
        super().plotImage(ax,image_names)

        ## colour bar
        if self.use_colorbar:
            addColorbar(
                ax,mcm.get_cmap(self.cmap),
                10**self.min_quantity,10**self.max_quantity,
                self.cbar_label,
                logflag = self.take_log_of_quantity,
                fontsize=self.fontsize,cmap_number=0)

def getImageGrid(
    BoxSize,
    Xmin,Xmax,
    Ymin,Ymax,
    Zmin,Zmax,
    npix_x,npix_y,
    pos,mass,quantity,
    take_log_of_quantity,
    hsml=None):

    ## set c-routine variables
    desngb   = 32
    Axis1    = 0
    Axis2    = 1
    Axis3    = 2

    Hmax     = 0.5*(Xmax-Xmin)

    n_smooth = pos.shape[0]

    ## output array for sum along the line of sight
    totalMassMap = np.zeros(shape = (npix_x,npix_y),dtype=np.float32)

    ## output array for average along the line of sight
    massWeightedQuantityMap = np.zeros(shape = (npix_x,npix_y),dtype=np.float32)
    
    ## create hsml output array
    if hsml is None:
        hsml = np.zeros(mass.shape[0],dtype=np.float32)
    else:
        print("Using provided smoothing lengths")
    
    c_f_p      = ctypes.POINTER(ctypes.c_float)
    pos_p      = pos.ctypes.data_as(c_f_p)
    hsml_p     = hsml.ctypes.data_as(c_f_p)
    mass_p     = mass.ctypes.data_as(c_f_p)
    quantity_p = quantity.ctypes.data_as(c_f_p)
    w_f_p    = totalMassMap.ctypes.data_as(c_f_p)
    q_f_p    = massWeightedQuantityMap.ctypes.data_as(c_f_p)

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

    # normalise by area of each pixel to get SFC density (column density)
    Acell = (Xmax-Xmin)/npix_x * (Ymax-Ymin)/npix_y
    columnDensityMap = totalMassMap/(Acell) # 10^10 Msun / kpc^-2 
    
    # convert into Msun/pc^2
    unitmass_in_g = 1.9890000e+43 
    solar_mass    = 1.9890000e+33
    conv_fac = (unitmass_in_g/solar_mass) / (1.0e3)**2 ## Msun/pc^2
    columnDensityMap *= conv_fac
    columnDensityMap = np.log10(columnDensityMap)
    print(
	'log10 minmax(columnDensityMap)',
	np.min(columnDensityMap),
	np.max(columnDensityMap))

    # massWeightedQuantityMap contains the mass-weighted quantity
    if take_log_of_quantity:
        massWeightedQuantityMap = np.log10(massWeightedQuantityMap)
    print(
	'log10 minmax(massWeightedQuantityMap)',
	np.min(massWeightedQuantityMap),
	np.min(massWeightedQuantityMap))
   
    return columnDensityMap,massWeightedQuantityMap



