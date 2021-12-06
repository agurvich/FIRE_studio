import os
import numpy as np 
import h5py

import matplotlib.pyplot as plt
from matplotlib.colors import rgb_to_hsv,hsv_to_rgb
#import matplotlib.gridspec as gridspec

from abg_python import append_function_docstring,filterDictionary
from abg_python.plot_utils import nameAxes

from abg_python.galaxy.gal_utils import Galaxy
from abg_python.galaxy.metadata_utils import Metadata,metadata_cache

from ..utils.stellar_utils.load_stellar_hsml import get_particle_hsml
from ..utils.camera_utils import Camera

class Drawer(object):

    def drawCoordinateAxes(
        self,
        ax,
        spacing=1,
        length=10,
        colors=None):

        if colors is None:
            colors = ['red','blue','green']
        
        ## create the axis lines
        points = np.arange(spacing,length+spacing,spacing)
        
        ## initialize the coordinate array, + 1 for the origin
        coordinates = np.zeros((3*points.size + 1,3))

        ## for each direction, fill one axis with the points
        for i in range(3):
            coordinates[points.size*i:points.size*(i+1),i] = points 
        
        ## perform the rotation
        coordinates = self.camera.rotate_array(coordinates,offset=True)

        ## plot the new x-y coordiantes
        for i in range(3):
            these_coords = coordinates[points.size*i:points.size*(i+1)]
            ax.plot(these_coords[:,0],these_coords[:,1],'.',color=colors[i])

        ax.plot(0,0,'.',c=colors[-1])
        for i,label in enumerate(['x','y','z']):
            x,y,z = coordinates[points.size*(i+1)-1]
            ax.text(x,y,label)

        return ax

    def plotImage(
        self,
        ax,
        final_image,
        cbar_label=None,
        **kwargs): 

        ## fill the pixels of the the scale bar with white
        if self.scale_bar:
            self.addScaleBar(final_image)

        ## main imshow call
        imgplot = ax.imshow(
            final_image, 
            extent = (self.Xmin,self.Xmax,self.Ymin,self.Ymax),
            origin = 'lower')

        ## set the y units to equal the x units-- make sure nothing
        ##  is being skewed. UNRELATED to self.aspect_ratio
        ax.set_aspect(1)

        ## turn off the axis if asked
        if self.noaxis:
            ax.axis('off')

        ## will check relevant flags internally
        self.addText(ax)

####### image utilities #######
    def addScaleBar(self,image):

        ## set scale bar length
        self.scale_label_text = r"$\mathbf{%1g \, \rm{kpc}}$"%self.scale_line_length

        # Convert to pixel space
        length_per_pixel = (self.Xmax - self.Xmin) / self.npix_x
        self.scale_line_length_px = int(self.scale_line_length / length_per_pixel)

        # Position in terms of image array indices
        scale_line_x_start = int(0.05 * self.npix_x)
        scale_line_x_end = min(scale_line_x_start + self.scale_line_length_px,self.npix_x)
        scale_line_y = int(0.02 * self.npix_y)

        npix_thick = 12
        # Go through pixels for scale bar, setting them to white
        for x_index in range(scale_line_x_start, scale_line_x_end):
            image[scale_line_y:scale_line_y+npix_thick, x_index,:3] = 1 if self.font_color in ['w','white'] else 0
        return image

    def addText(self,ax):
        ## handle any text additions
        if self.figure_label is not None:
        ## plot the  figure label in the top right corner
            nameAxes(ax,None,None,None,
                supertitle=self.figure_label,
                subfontsize=self.fontsize,
                font_color=self.font_color,
                swap_annotate_side=self.figure_label_side=='right')

        ## Set parameters

        scale_label_position = 0.06 

        if self.scale_bar: 
            ## plot the scale bar label
            label2 = ax.text(
                scale_label_position,0.03,
                self.scale_label_text,
                fontweight = 'bold',
                fontsize=self.fontsize*0.75,
                transform = ax.transAxes,
                verticalalignment='bottom')
            label2.set_color(self.font_color)

    def renormalizeTransposeImage(self,image,min_val,max_val,quantity_name):
        if self.master_loud:
            print('min_%s = '%quantity_name,min_val)
            print('max_%s = '%quantity_name,max_val)

            print('Image range (%s): '%quantity_name,np.nanmin(image),np.nanmax(image))
        image = image - min_val
        image = image / (max_val - min_val)
        
        ## clip anything outside the range
        image[image < 0.0] = 0.0
        image[image > 1.0] = 1.0
        image = image*255.0

        if self.master_loud:
            print('Image range (8bit): ',np.nanmin(image),np.nanmax(image))

        ## cast to integer to use as indices for cmap array
        image = image.astype(np.uint16) 

        return image.T

    def saveFigure(
        self,
        ax,
        image_name=None,
        **savefig_args):

        if self.noaxis:
            ## remove whitespace around the axis, apparently the x/y origin is offset in pixel 
            ## space and so the bounding box doesn't actually reflect the left/bottom edge of the 
            ## axis
            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.yaxis.set_major_locator(plt.NullLocator())
            savefig_args['bbox_inches']='tight'
            savefig_args['pad_inches']=0

        
        if image_name is None:
            image_name = "%03d_%dkpc.pdf" % (self.snapnum, 2*self.camera.camera_dist)

        if 'png' not in image_name and 'pdf' not in image_name:
            image_name+='.pdf'

        ax.get_figure().savefig(
            os.path.join(self.datadir,image_name),
            dpi=300,
            **savefig_args)
    
    def gradientBlendImages(self,image_1,image_2=None,**kwargs):

        if image_2 is None: image_2 = self.produceImage(**kwargs)

        new_image = np.zeros(image_1.shape)

        xs = np.ones(new_image.shape[0])*0.05

        gradient = np.linspace(0.05,1,xs.shape[0]//2)[::-1]
        xs[:xs.shape[0]//2] = gradient


        new_image = image_1*xs[:,None]+image_2*(1-xs[:,None])
        hsv_image = rgb_to_hsv(new_image)
        hsv_image[...,-1]*=1.2
        new_image = hsv_to_rgb(hsv_image)

        return new_image

class Studio(Drawer):
    """`FIREstudio` parent class that regularizes image setup, rotation, 
        caching, etc between `GasStudio` and `StarStudio` classes. 
    """

    def __init__(
        self,
        datadir, ## directory to put intermediate and output files 
        snapnum, ## snapshot number 
        sim_name,
        cache_file_name=None, ##  the name of the file to save maps to
        gas_snapdict=None, ## open galaxy gas snapshot dictionary 
        star_snapdict=None, ## open galaxy star snapshot dictionary
        galaxy_kwargs=None,
        master_loud=True,
        **kwargs
        ):
        """Initializes a cache file to be read from and sets any image parameters
            that have been passed to init. Extra `kwargs` will be passed to `set_ImageParams`.
            Snapshot dictionaries can be passed to avoid `FIREstudio` using `load_SnapshotData`
            to open them from disk. They require a minimum number of keys to function.

`GasStudio` requires: 
```python
    snapdict['Coordinates'] ## coordinates of the particles
```
(and ideally `'SmoothingLengths'`, in the same units as coordinates, but these can be calculated). 

`StarStudio` requires: 
```python
gas_snapdict['Coordinates'] ## coordinates of the particles
gas_snapdict['Metallicity'] ## metallicity (mass fractions) of the particles 
gas_snapdict['Masses'] ## masses of the particles in 1e10 solar masses
gas_snapdict['Temperature'] ## temperature of the gas in K

star_snapdict['Coordinates'] ## coordinates of the particles
star_snapdict['Metallicity'] ## metallicity (mass fractions) of the particles 
star_snapdict['Masses'] ## masses of the particles in 1e10 solar masses
star_snapdict['AgeGyr'] ## age of particles in Gyr
```
(and ideally `'SmoothingLengths'` for both, in the same units as coordinates for both, but these can be calculated). 

            Input:

                datadir -- directory to put intermediate and output files 
                snapnum -- snapshot number (feel free to lie if you aren't using
                    FIRE_studio to open a snapshot, it is needed for cache file name though)
                sim_name -- name of the simulation, i.e. m12i_res7100. prepends the cache_file_name
                    if the sim_name isn't already in the path to disambiguate caches.  
                cache_file_name = None -- the name of the file to save maps to
                    defaults to proj_maps_%03d.hdf5, changes when you use set_DataParams
                    to choose a snapshot.

                gas_snapdict = None -- an open galaxy gas snapshot dictionary, opens from snapdir/snapnum if None
                    (must then also provide snapdir)
                star_snapdict = None -- an open galaxy star snapshot dictionary, opens from snapdir/snapnum if None 
                    (must then also provide snapdir)
                galaxy_kwargs = None -- dictionary that contains kwargs that should be passed to the opened
                    abg_python.galaxy.Galaxy instance that is used to load snapshot data from disk. 
                master_loud = True -- flag for turning off *all* print statements

                **set_ImageParams_kwargs -- keyword arguments that will be passed to
                    set_ImageParams

            Output:

                None"""
        
        ## bind loud flag
        self.master_loud = master_loud

        if cache_file_name is None:
            cache_file_name = "proj_maps"

        self.cache_file_name = cache_file_name
        self.cache_file = None
        
        ## IO stuff
        ##  append 'firestudio' if we're living in a folder that belongs to a Galaxy
        if 'firestudio' not in datadir and sim_name in datadir:
            datadir = os.path.join(datadir,'firestudio')
        self.datadir = datadir

        ## make the datadirectory if it doesn't exist
        if not os.path.isdir(self.datadir):
            os.makedirs(self.datadir)

        self.gas_snapdict = gas_snapdict
        self.star_snapdict = star_snapdict
        
        ## {}.update enforces that galaxy_kwargs is a dictionary, lol...
        self.galaxy_kwargs = {} if galaxy_kwargs is None else {}.update(galaxy_kwargs)

        ## create, if necessary, directories to store intermediate and output files,
        ##  this could get crowded! sets self.image_dir and self.projection_dir
        #self.makeOutputDirectories(datadir)
 
        if 'camera' not in kwargs or kwargs['camera'] is None:
            camera_kwargs = {'camera_pos':[0,0,15],'camera_focus':[0,0,0],'camera_north':None}
            for ckwarg in list(camera_kwargs.keys())+['quaternion']:
                if ckwarg in kwargs: camera_kwargs[ckwarg] = kwargs.pop(ckwarg)
            kwargs['camera'] = Camera(**camera_kwargs)

        ## initialize the object with some default image params that will
        ##  make a face-on image of a galaxy, can always set these manually
        ##  with  external calls to set_ImageParams
        self.set_ImageParams(use_defaults=True,snapnum=snapnum,sim_name=sim_name,**kwargs)

####### I/O functions #######
    def load_SnapshotData(
        self,
        gas_mask=None,
        star_mask=None,
        **kwargs):
        """Binds simulation output to self.gas_snapdict and self.star_snapdict.

            Input:

                use_saved_subsnapshots = False -- save/load subsnapshots, uncompressed copies of the snapshot
                    oriented on the main disk with particles within the virial radius. This can 
                    take up lots of disk space. 
                del_galaxy = True -- whether the abg_python.galaxy.Galaxy object should be deleted after 
                    being used to get the snapshot dictionaries.

            Output:

                None/abg_python.galaxy.Galaxy object -- if del_galaxy == False then returns the galaxy object, 
                    otherwise returns None. """

        ## determine if we need to open any snapshot data
        if (self.gas_snapdict is None or 
            self.gas_snapdict['snapnum'] != self.snapnum ): ## haven't loaded this data yet, or we are replacing it
            return_value = self.__get_snapdicts(
                self.sim_name,
                self.snapdir,self.snapnum,**kwargs)

            if hasattr(self,'masked_gas_snapdict'):
                del self.masked_gas_snapdict

            if hasattr(self,'masked_star_snapdict'):
                del self.masked_star_snapdict
    

        ## apply a mask to the snapdict if requested to
        if gas_mask is not None:
            self.masked_gas_snapdict = filterDictionary(self.gas_snapdict,gas_mask)
        if star_mask is not None:
            self.masked_star_snapdict = filterDictionary(self.star_snapdict,star_mask)

        return return_value

    def get_HSML(
        self,
        snapdict_name,
        use_metadata=True,
        save_meta=True,
        assert_cached=False,
        loud=True,
        **kwargs, 
        ):
        """Compute smoothing lengths for particles that don't have them,
            typically collisionless particles (like stars). 

            Input:

                snapdict_name -- name in the form of `'%s_snapdict'%snapdict_name`
                    that will be used to compute smoothing lengths for. 

                use_metadata = True -- flag to search cache for result
                save_meta = True -- flag to cache the result
                assert_cached = False -- flag to require a cache hit
                loud = True -- whether cache hits/misses should be announced
                    to the console.
                
            Output:

                smoothing_lengths -- numpy array of estimated smoothing lengths"""

        @metadata_cache(
            '%s_data'%snapdict_name,  ## hdf5 file group name
            ['%s_SmoothingLengths'%snapdict_name],
            use_metadata=use_metadata,
            save_meta=save_meta,
            assert_cached=assert_cached,
            loud=loud,
            force_from_file=True) ## read from cache file, not attribute of object
        def compute_HSML(self):
            snapdict = getattr(self,snapdict_name+'_snapdict')
            pos = snapdict['Coordinates']
            smoothing_lengths = get_particle_hsml(pos[:,0],pos[:,1],pos[:,2])
            return smoothing_lengths

        return compute_HSML(self,**kwargs)

    def __get_snapdicts(
        self,
        sim_name,
        snapdir,snapnum,
        use_saved_subsnapshots=False,
        del_galaxy=True,
        **kwargs):

        these_kwargs = self.galaxy_kwargs.copy()
        these_kwargs.update(kwargs)

        galaxy = Galaxy(
            sim_name,
            snapdir,
            snapnum,
            datadir=os.path.dirname(self.datadir), 
            loud_metadata=False, ## shh don't let them know
            save_header_to_table=False,
            **these_kwargs)

        ## add use_saved_subsnapshots into these kwargs to be passed to extract
        ##  main halo to allow users to use cached snapshot data
        these_kwargs.update({'use_saved_subsnapshots':use_saved_subsnapshots})

        ## handles opening the snapshot, centering it, and rotating it to be face-on.
        galaxy.extractMainHalo(
            save_meta=False,
            **these_kwargs) ## metadata cache will pull only the good keys out

        ## bind the snapshot dictionaries
        self.gas_snapdict = galaxy.sub_snap
        self.star_snapdict = galaxy.sub_star_snap

        ## just get rid of it now that we've opened it
        if del_galaxy:
            del galaxy
        else:
            return galaxy

    def set_ImageParams(
        self,
        this_setup_id=None,
        use_defaults=False,
        loud=True,
        **kwargs):
        """Changes the parameters of the image. If `use_defaults=True` then 
            default values of the parameters will be set. Leave `use_defaults=False`
            to adjust only the keywords passed. 

            Input:

                this_setup_id = None -- string that identies this image setup in the cache,
                    use this to differentiate between things that might otherwise
                    be overwritten.
                use_defaults = False -- flag to overwrite any kwargs that *aren't* passed in
                    this call to their default value. 
                loud = True -- 

                frame_half_thickness = 15 -- half-thickness of image in z direction

                aspect_ratio = 1 -- shape of image, y/x TODO figure out if this is necessary to pass?
                pixels = 1200 -- pixels in x direction, resolution of image
                figure_label = '' -- string to be put in upper right corner
                'figure_label_side' = 'right' --  corner to put label in

                scale_bar = True -- flag to plot length scale bar in lower left corner
                scale_line_length = 5 -- length of the scale bar in kpc

                noaxis = True -- turns off axis ticks
                savefig = None -- save the image as a png if passed a string
                fontsize = 12  -- fontsize of figure label and scale bar text
                font_color = 'white' -- color of the subtitle font
                
                snapdir - path to simulation output
                snapnum -  which snapshot to open
                sim_name - name of simulation (i.e. m12i_res7100)

            Output:

                None

Example usage:
```python
studio.set_ImageParams(
    this_setup_id='my_custom_setup',
    scale_bar=False,
    figure_label='high redshift')
```"""

        default_kwargs = {
            'camera':None,
            'frame_half_thickness':None, ## half-thickness of image in z direction
            'aspect_ratio':1, ## shape of image, y/x TODO figure out if this is necessary to pass?
            'pixels':1200, ## pixels in x direction, resolution of image
            'figure_label':'', ## string to be put in upper right corner
            'figure_label_side':'right', ## corner to put label in
            'scale_bar':True,  ## flag to plot length scale bar in lower left corner
            'scale_line_length':5, ## length of the scale line in kpc
            'noaxis':True, ## turns off axis ticks
            'savefig':None, ## save the image as a png if passed a string
            'fontsize':12,  ## font size of scale bar and figure label
            'font_color':'w', ## font color of scale bar and figure label
            'snapdir':None,
            'snapnum':None,
            'sim_name':None
            }

        for kwarg in kwargs:
            ## only set it here if it was passed
            if kwarg in default_kwargs:
                value = kwargs[kwarg]
                ## remove it from default_kwargs
                default_kwargs.pop(kwarg)

                ## set it to the object
                if loud and self.master_loud:
                    print("setting",kwarg,
                        'to user value of:',value)
                ## set it to the object
                setattr(self,kwarg,value)
            else:
                if self.master_loud:
                    print(kwarg,'ignored. Did you mean something else?',
                        default_kwargs.keys())

        if use_defaults:
            ## if we haven't already set the camera by passed kwarg
            ##  need to replace None with [0,0,0]
            if ('camera' in default_kwargs and 
                default_kwargs['camera'] is None):
                raise ValueError('Cannot explicitly set camera = None,'+
                ' pass a utils.camera_utils.Camera instance instead.')

            if ('frame_half_thickness' in default_kwargs and 
                default_kwargs['frame_half_thickness'] is None):
                default_kwargs['frame_half_thickness'] = self.camera.camera_dist
             
            ## set the remaining image parameters to their default values
            for default_arg in default_kwargs:
                value = default_kwargs[default_arg]
                if loud and self.master_loud:
                    print("setting",default_arg,
                        'to default value of:',value)
                setattr(self,default_arg,value)

        ## determine the edges of our frame so we can cull the rest later,
        ##  also set cell area info
        self.computeFrameBoundaries()
        
        ## identify the combination of these parameters
        if this_setup_id is None:
            this_setup_id = self.__identifyThisSetup()

        self.this_setup_id = this_setup_id

        self.set_CacheFile()

    set_ImageParams.default_kwargs = {
            'camera':None,
            'frame_half_thickness':None, ## half-thickness of image in z direction
            'aspect_ratio':1, ## shape of image, y/x TODO figure out if this is necessary to pass?
            'pixels':1200, ## pixels in x direction, resolution of image
            'figure_label':'', ## string to be put in upper right/left corner
            'figure_label_side':'right', ## corner to put label in
            'scale_bar':True,  ## flag to plot length scale bar in lower left corner
            'scale_line_length':5, ## length of the scale line in kpc
            'noaxis':True, ## turns off axis ticks
            'savefig':None, ## save the image as a png if passed a string
            'fontsize':12,  ## font size of scale bar and figure label
            'font_color':'w',
            'snapdir':None,
            'snapnum':None,
            'sim_name':None}

    def set_CacheFile(self):
        """ Creates the cache hdf5 file. Requires self.snapnum and sim_name be set.

            Input:

                None

            Output:

                cache_file -- abg_python.galaxy.metadata_utils.Metadata object"""

        if (self.snapnum is not None and 
            self.sim_name is not None):
            ## read the snapshot number from the end of the metapath
            if (self.cache_file is None or 
                int(self.cache_file.metapath.split('_')[-1].split('.hdf5')[0]) != self.snapnum):
                h5name=self.cache_file_name+"_%03d.hdf5"%self.snapnum
                
                ## disambiguate to ensure simulation data never gets mixed
                if (self.sim_name not in self.datadir and
                    self.sim_name not in h5name):
                    h5name = self.sim_name+'_'+h5name

                self.metadata = self.cache_file = Metadata(os.path.join(self.datadir,h5name)) 
        else:
            raise IOError("Need to set self.snapnum and self.sim_name to disambiguate cache_file")

        return self.metadata

    def print_ImageParams(self):
        """ Prints current image setup to console.

            Input:

                None

            Output:

                None"""
        default_kwargs = [
            'frame_half_thickness',
            'aspect_ratio',
            'pixels', 
            'figure_label', 
            'figure_label_side',
            'scale_bar',  
            'noaxis',
            'savefig',
            'fontsize',
            'snapdir',
            'snapnum']

        ## print the current value, not the default value
        for arg in default_kwargs:
            print(arg,'=',getattr(self,arg))

    def __identifyThisSetup(self):
        ## uniquely identify this projection setup using a simple "hash"
        rotation_string = 'quat_%.2f_%.2f_%.2f_%.2f'%(
            np.round(self.camera.quaternion[0],decimals=2),
            np.round(self.camera.quaternion[1],decimals=2),
            np.round(self.camera.quaternion[2],decimals=2),
            np.round(self.camera.quaternion[3],decimals=2))

        self.this_setup_id = (
        "npix%d_width%.2fkpc_depth%.2fkpc_x%.2f_y%.2f_z%.2f_%s_aspect%.2f"%(
            self.pixels, 
                np.round(self.camera.camera_dist,decimals=2),
                np.round(self.frame_half_thickness,decimals=2),
                np.round(self.camera.camera_focus[0],decimals=2),
                np.round(self.camera.camera_focus[1],decimals=2),
                np.round(self.camera.camera_focus[2],decimals=2),
                rotation_string,
                np.round(self.aspect_ratio,decimals=2)))
        return self.this_setup_id

    def computeFrameBoundaries(self):
        """Uses current image parameters to calculate the minimum and maximum 
            x, y, and z limits. 

            Input: 
                
                None

            Output:

                None

            Attributes set:
                self.Xmin, self.Xmax -- 
                self.Ymin, self.Ymax -- 
                self.Zmin, self.Zmax -- 
                self.npix_x, self.npix_y -- 
                self.Acell -- 
                
        """
        ## +- camera_dist limits -> 45 degree FOV
        self.Xmin,self.Xmax = self.camera.camera_focus[0] + np.array(
            [-self.camera.camera_dist,self.camera.camera_dist])

        self.Ymin,self.Ymax = self.camera.camera_focus[1] + np.array(
            [-self.camera.camera_dist,self.camera.camera_dist])*self.aspect_ratio

        self.Zmin,self.Zmax = self.camera.camera_focus[2] + np.array(
            [-self.frame_half_thickness,self.frame_half_thickness])

        ## Set image size 
        self.npix_x   = self.pixels #1200 by default
        self.npix_y   = int(self.pixels*self.aspect_ratio) #1200 by default

        self.Acell = (self.Xmax-self.Xmin)/self.npix_x * (self.Ymax-self.Ymin)/self.npix_y

    def cullFrameIndices(self,Coordinates):

        ## extract a cube of particles that are in relevant area
        ind_box = ((Coordinates[:,0] > self.Xmin) & (Coordinates[:,0] < self.Xmax) &
                   (Coordinates[:,1] > self.Ymin) & (Coordinates[:,1] < self.Ymax) &
                   (Coordinates[:,2] > self.Zmin) & (Coordinates[:,2] < self.Zmax))

        return ind_box

## append method docstrings to class docstring
append_function_docstring(Studio,Studio.__init__)
append_function_docstring(Studio,Studio.load_SnapshotData)
append_function_docstring(Studio,Studio.print_ImageParams)
append_function_docstring(Studio,Studio.computeFrameBoundaries)
append_function_docstring(Studio,Studio.get_HSML)
