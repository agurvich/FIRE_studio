import os
import numpy as np 

import matplotlib.pyplot as plt
from matplotlib.colors import rgb_to_hsv,hsv_to_rgb
#import matplotlib.gridspec as gridspec

from abg_python import append_function_docstring,filterDictionary
from abg_python.plot_utils import nameAxes

from abg_python.galaxy.gal_utils import Galaxy
from abg_python.galaxy.metadata_utils import Metadata,metadata_cache

from ..utils.stellar_utils.load_stellar_hsml import get_particle_hsml
from ..utils.camera_utils import Camera

try:
    from numba import njit
    @njit()
    def sample_gradient_at_angle(gradient_mask,image_1,angle):
        new_gradient_mask = np.zeros(image_1.shape[:-1])
        ## coordinates are w.r.t. center of image so rotation is correct
        for pix_i in range(-image_1.shape[0]//2,image_1.shape[0]//2):
            for pix_j in range(-image_1.shape[1]//2,image_1.shape[1]//2):
                lookup_i = int(np.round(np.cos(angle)*pix_i-np.sin(angle)*pix_j,0))+gradient_mask.shape[0]//2
                lookup_j = int(np.round(np.sin(angle)*pix_i+np.cos(angle)*pix_j,0))+gradient_mask.shape[1]//2
                new_gradient_mask[pix_i+image_1.shape[0]//2,pix_j+image_1.shape[1]//2] = gradient_mask[lookup_i,lookup_j]
        return new_gradient_mask
except ImportError:
    ## don't have numba so gradient blending will be slow. oh well
    def sample_gradient_at_angle(gradient_mask,image_1,angle):
        new_gradient_mask = np.zeros(image_1.shape[:-1])
        ## coordinates are w.r.t. center of image so rotation is correct
        for pix_i in range(-image_1.shape[0]//2,image_1.shape[0]//2):
            for pix_j in range(-image_1.shape[1]//2,image_1.shape[1]//2):
                lookup_i = int(np.round(np.cos(angle)*pix_i-np.sin(angle)*pix_j,0))+gradient_mask.shape[0]//2
                lookup_j = int(np.round(np.sin(angle)*pix_i+np.cos(angle)*pix_j,0))+gradient_mask.shape[1]//2
                new_gradient_mask[pix_i+image_1.shape[0]//2,pix_j+image_1.shape[1]//2] = gradient_mask[lookup_i,lookup_j]
        return new_gradient_mask


class Drawer(object):
    def render(
        self,
        ax:plt.Axes=None,
        **kwargs):
        """ Generates an image with the `produceImage` method and then plots it with the `plotImage` method.

            Input: 

                ax = None -- axis to plot image to, if None will create a new figure

            Output:

                ax -- the axis the image was plotted to
                final_image -- Npixels x Npixels x 3 RGB pixel array"""

        if ax is None:
            fig,ax = plt.figure(),plt.gca()
        else:
            fig = ax.get_figure()

        ## remap the C output to RGB space
        final_image = self.produceImage(**kwargs)

        ## plot that RGB image and overlay scale bars/text
        self.plotImage(ax,final_image)

        ## save the image
        if self.savefig is not None:
            self.saveFigure(fig,self.savefig)

        return ax,final_image

    def drawCoordinateAxes(
        self,
        ax:plt.Axes,
        spacing:float=1,
        length:float=10,
        colors:list=None):
        """[summary]

        Parameters
        ----------
        ax : plt.Axes
            [description]
        spacing : float, optional
            [description], by default 1
        length : float, optional
            [description], by default 10
        colors : list, optional
            [description], by default None

        Returns
        -------
        [type]
            [description]
        """

        if colors is None:
            colors = ['red','green','blue']
        
        ## create the axis lines
        points = np.arange(spacing,length+spacing,spacing)
        
        ## initialize the coordinate array, + 1 for the origin
        coordinates = np.zeros((3*points.size + 1,3))

        ## for each direction, fill one axis with the points
        for i in range(3):
            coordinates[points.size*i:points.size*(i+1),i] = points 
        
        ## perform the rotation
        coordinates,mask = self.camera.project_and_clip(coordinates)
        if np.sum(mask) == 0: return ax

        ## plot the new x-y coordiantes
        for i in range(3):
            these_coords = coordinates[points.size*i:points.size*(i+1)]
            ax.plot(these_coords[:,0],these_coords[:,1],'.',color=colors[i])

        ax.plot(0,0,'.',c=colors[-1])
        for i,label in enumerate(['x','y','z']):
            x,y,z = coordinates[points.size*(i+1)-1]
            ax.text(x,y,label,fontdict={'color':'white'})

        ax.set_facecolor('k')
        ax.set_aspect(1)
        ax.set_xlim(-length,length)
        ax.set_ylim(-length,length)

        subtitle = str(np.round(self.camera.project_array(np.identity(3),offset=False),1))
        subtitle = subtitle.replace('[','').replace(']','')
        nameAxes(
            ax,None,None,None,
            subtitle=subtitle,
            supertitle=str(np.round(self.camera.camera_pos,0)),
            swap_annotate_side=True,font_color='w')
        return ax

    def plotImage(
        self,
        ax:plt.Axes,
        final_image:np.ndarray,
        **kwargs): 
        """Bsae method for overlaying artists on top of projected image.
        if `self.scale_bar`:
            overlays a scale bar by filling the RGB pixel values with white
        if `self.noaxis`:
            removes the coordinate axes, labels, and ticks

        will also add `self.figure_label` as text to the image. 
        See `~firestudio.studios.studio.Studio.set_ImageParams` for details.

        Parameters
        ----------
        ax : plt.Axes
            axis to plot image to 
        final_image : np.ndarray
            array of RGB image pixel values
        """

        ## fill the pixels of the the scale bar with white
        if self.scale_bar:
            self.addScaleBar(final_image)

        ## main imshow call
        imgplot = ax.imshow(
            final_image, 
            extent = (
                -self.camera.frame_half_width,
                self.camera.frame_half_width,
                -self.camera.frame_half_width,
                self.camera.frame_half_width),
            origin = 'lower')

        ## set the y units to equal the x units-- make sure nothing
        ##  is being skewed. UNRELATED to self.aspect_ratio
        ax.set_aspect(1)

        ## turn off the axis if asked
        if self.noaxis:
            ax.axis('off')
            ## remove whitespace around the axis, apparently the x/y origin is offset in pixel 
            ## space and so the bounding box doesn't actually reflect the left/bottom edge of the 
            ## axis
            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.yaxis.set_major_locator(plt.NullLocator())

        ## will check relevant flags internally
        self.addText(ax)

####### image utilities #######
    def addScaleBar(self,image:np.ndarray):
        """[summary]

        Parameters
        ----------
        image : np.ndarray
            [description]

        Returns
        -------
        [type]
            [description]
        """

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

    def addText(self,ax:plt.Axes):
        """[summary]

        Parameters
        ----------
        ax : plt.Axes
            [description]
        """
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

    def renormalizeTransposeImage(
        self,
        image:np.ndarray,
        min_val:float,
        max_val:float,
        quantity_name:str):
        """[summary]

        Parameters
        ----------
        image : np.ndarray
            [description]
        min_val : float
            [description]
        max_val : float
            [description]
        quantity_name : str
            [description]

        Returns
        -------
        [type]
            [description]
        """

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
        fig,
        image_name:str=None,
        **savefig_args):
        """

        Parameters
        ----------
        fig : [type]
            [description]
        image_name : str, optional
            [description], by default None
        """

        if self.noaxis:
            savefig_args['bbox_inches']='tight'
            savefig_args['pad_inches']=0
        
        if image_name is None:
            image_name = "%03d_%dkpc.pdf" % (self.snapnum, 2*self.camera.camera_dist)

        if 'png' not in image_name and 'pdf' not in image_name:
            image_name+='.pdf'

        fig.savefig(
            os.path.join(self.datadir,image_name),
            dpi=300,
            **savefig_args)
    
    def gradientBlendImages(
        self,
        image_1:np.ndarray,
        image_2:np.ndarray=None,
        gradient_width_percent:float=0.1,
        angle:float=None,
        **kwargs):
        """[summary]

        Parameters
        ----------
        image_1 : np.ndarray
            [description]
        image_2 : np.ndarray, optional
            [description], by default None
        gradient_width_percent : float, optional
            [description], by default 0.1
        angle : float, optional
            [description], by default None

        Returns
        -------
        [type]
            [description]
        """

        if image_2 is None: image_2 = self.produceImage(**kwargs)

        ## crop images as necessary
        image_1 = image_1[:image_2.shape[0],:image_2.shape[1]]
        image_2 = image_2[:image_1.shape[0],:image_1.shape[1]]

        if image_1.shape[-1] == 3:
            print('appending')
            image_1 = np.append(image_1,np.ones(image_1.shape[:-1])[...,None],axis=-1)
        if image_2.shape[-1] == 3:
            print('appending')
            image_2 = np.append(image_2,np.ones(image_2.shape[:-1])[...,None],axis=-1)

        ## build the gradient to sample from
        gradient_edge_length = np.ceil(np.sqrt(2*image_1.shape[0]*image_1.shape[1])).astype(int)
        gradient_mask = np.zeros((gradient_edge_length,gradient_edge_length))
        gradient_mask[:gradient_mask.shape[0]//2,:] = 1
        
        ## TODO could consider allowing gradient to be offcenter I suppose
        offset = (gradient_mask.shape[0]//2-image_1.shape[0]//2)
        gradient_begin_index = int(image_1.shape[0]*(0.5-gradient_width_percent/2))+offset
        gradient_end_index = int(image_1.shape[0]*(0.5+gradient_width_percent/2))+offset
        gradient_mask[gradient_begin_index:gradient_end_index] = np.linspace(1,0,gradient_end_index-gradient_begin_index)[:,None]
        
        ## expect angle in degrees
        if angle is not None: angle*=np.pi/180
        
        ## sample gradient but rotate coordinates. this takes a while because it's a nested for loop over pixels
        ##  so we can use numba if we have it
        new_gradient_mask = sample_gradient_at_angle(gradient_mask,image_1,angle)
        
        final_image = np.ones(image_1.shape)
        ## apply the alpha blend from the gradient
        final_image[...,:-1] = image_1[...,:-1]*(new_gradient_mask*image_1[...,-1])[...,None]+image_2[...,:-1]*(((1-new_gradient_mask)*image_2[...,-1])[...,None])

        ## final alpha channel is 1s everywhere
        return final_image
class Studio(Drawer):
    """
    `FIREstudio` parent class that regularizes image setup, rotation, 
        caching, etc between `GasStudio` and `StarStudio` classes. 

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
    (and ideally `'SmoothingLengths'` for both, in the same units as coordinates, but these can be calculated)."""

    def __repr__(self):
        """ implementation of built-in __repr__ method for printing.

        Returns
        -------
        str
            string to print back to the console
        """
        return 'Studio instance'

    def __init__(
        self,
        datadir:str,  
        snapnum:int,  
        sim_name:str,
        cache_file_name:str=None, 
        gas_snapdict:dict=None,  
        star_snapdict:dict=None, 
        galaxy_kwargs:dict=None,
        master_loud:bool=True,
        setup_id_append:str='',
        **kwargs
        ): 
        """ Base class that handles camera manipulation and data caching. 

        Parameters
        ----------
        datadir : str
            directory to put intermediate and output files, 'firestudio' is appended
            if the directory contains sim_name
        snapnum : int
            snapshot number (feel free to lie if you aren't using
            FIRE_studio to open a snapshot, it is needed for cache file name though)
        sim_name : str
            name of the simulation, i.e. m12i_res7100. prepends the cache_file_name
            if the sim_name isn't already in the path to disambiguate caches.  
        cache_file_name : str, optional
            the name of the file to save maps to,
            by default 'proj_maps_%03d.hdf5'%snapnum
        gas_snapdict : dict, optional
            a dictionary containing gas data (or that which should be treated as 
            gas data, depending on the context), by default None
        star_snapdict : dict, optional
            a dictionary containing collisionless particle data (i.e. star particles, or
            that which should be treated as star data, depending on the context),
            by default None
        galaxy_kwargs : dict, optional
            dictionary that contains kwargs that should be passed to the opened
            abg_python.galaxy.Galaxy instance that is used to load snapshot data from disk,
            by default None
        master_loud : bool, optional
            flag for enabling/disabling *all* print statements, by default True
        """
        
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
            camera_kwargs = {
                'camera_pos':[0,0,15],
                'camera_focus':[0,0,0],
                'camera_north':None}
            for ckwarg in list(camera_kwargs.keys())+['quaternion']:
                if ckwarg in kwargs: camera_kwargs[ckwarg] = kwargs.pop(ckwarg)
            kwargs['camera'] = Camera(**camera_kwargs)

        ## initialize the object with some default image params that will
        ##  make a face-on image of a galaxy, can always set these manually
        ##  with  external calls to set_ImageParams
        self.set_ImageParams(use_defaults=True,snapnum=snapnum,sim_name=sim_name,**kwargs)
        self.this_setup_id+=setup_id_append

####### I/O functions #######
    def load_SnapshotData(
        self,
        gas_mask:np.ndarray=None,
        star_mask:np.ndarray=None,
        **kwargs):
        """ Binds simulation output to self.gas_snapdict and self.star_snapdict.

        Parameters
        ----------
        gas_mask : bool np.ndarray, optional
            boolean mask that should be applied to the galaxy.sub_snap, by default None
        star_mask : bool np.ndarry, optional
            boolean mask that should be applied to the galaxy.sub_star_snap, by default None
        
        Keywords
        --------
            use_saved_subsnapshots: bool, optional
                save/load subsnapshots, uncompressed copies of the snapshot
                oriented on the main disk with particles within the virial radius. This can 
                take up lots of disk space, by default False
            del_galaxy: bool, optional 
                flag for whether the abg_python.galaxy.gal_utils.Galaxy object should be deleted after 
                being used to get the snapshot dictionaries.

        Returns
        -------
        ``abg_python.galaxy.Galaxy``
            None/abg_python.galaxy.Galaxy object -- if del_galaxy == False then returns the galaxy object, 
            otherwise returns None.

        """

        

        ## determine if we need to open any snapshot data
        if (self.gas_snapdict is None or 
            self.gas_snapdict['snapnum'] != self.snapnum ): ## haven't loaded this data yet, or we are replacing it
            return_value = self.__get_snapdicts(**kwargs)

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
        snapdict_name:str,
        use_metadata:bool=True,
        save_meta:bool=True,
        assert_cached:bool=False,
        loud:bool=True,
        **kwargs, 
        ):
        """ Compute smoothing lengths for particles that don't have them,
            typically collisionless particles (like stars). 

        Parameters
        ----------
        snapdict_name : str
            string identifying which snapdict should be used
            used to compute smoothing lengths, either 'gas' or 'star'
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
        np.float32 np.ndarray
            estimated smoothing lengths
        """

        

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
        use_saved_subsnapshots:bool=False,
        del_galaxy:bool=True,
        **kwargs):
        """ Open an abg_python.galaxy.gal_utils.Galaxy instance to load
            snapshot data from disk.

        Parameters
        ----------
        use_saved_subsnapshots : bool, optional
            flag to cache particle data within the virial radius by 
            saving it to disk in a "subsnapshot" hdf5 file, by default False
        del_galaxy : bool, optional
            delete the Galaxy object after binding the gas and star snapdicts,
            else return it, by default True

        Returns
        -------
        abg_python.galaxy.gal_utils.Galaxy
            the Galaxy instance used to open snapshot data, 
            only returned if `del_galaxy=False`
        """

        these_kwargs = self.galaxy_kwargs.copy()
        these_kwargs.update(kwargs)

        galaxy = Galaxy(
            self.sim_name,
            self.snapnum,
            datadir=os.path.dirname(self.datadir), 
            loud_metadata=False, ## shh don't let them know
            save_header_to_table=False,
            **these_kwargs)

        ## handles opening the snapshot, centering it, and rotating it to be face-on.
        galaxy.extractMainHalo(
            save_meta=False,
            use_saved_subsnapshots=use_saved_subsnapshots,
            **these_kwargs) ## metadata cache will pull only the good keys out

        ## bind the snapshot dictionaries
        self.gas_snapdict = galaxy.sub_snap
        self.star_snapdict = galaxy.sub_star_snap

        ## just get rid of it now that we've opened it
        if del_galaxy: del galaxy
        else: return galaxy

    def set_ImageParams(
        self,
        this_setup_id:str=None,
        use_defaults:bool=False,
        loud:bool=True,
        **kwargs):
        """ Changes the parameters of the image such as camera orientation, frame size, etc. 
            If `use_defaults=True` then default values of the parameters will be set and will 
            overwrite the current state. Leave `use_defaults=False` to adjust only the keywords passed.

        Parameters
        ----------
        this_setup_id : str, optional
            string to use to identify this combination of image parameters. If None, 
            then the image parameters are stringified and combined, by default None
        use_defaults : bool, optional
            overwrite current state with default values of each parameter, useful for initialization
            or resetting after making changes, by default False
        loud : bool, optional
            flag to print which parameters are being set/updated, by default True
        
        Keywords
        --------
        aspect_ratio: float
            ratio of number of pixels in each direction determining 
            the shape of image, y/x, by default 1
        pixels: int
            pixels in x direction, resolution of image, by default 1200
        figure_label: str
            string to be put in upper right corner, by default ''
        figure_label_side: str 
            side of the image to put label in, by default 'right' 
        scale_bar: bool
            flag to plot length scale bar in lower left corner, by default True
        scale_line_length: float
            length of the scale bar in kpc, by default 5
        noaxis: bool
            turns off axis ticks and labels, by default True
        savefig: str
            save the image as a png if passed a string or does not save a figure if None, by default None
        fontsize: int
            fontsize (in pt) of figure label and scale bar text, 12
        font_color: str/RGBA tuple
            color of the subtitle font, by default 'white'
        snapdir: str
            path to simulation output
        snapnum: int
            which snapshot to open/use for naming the cache
        sim_name: str
            name of simulation (i.e. m12i_res7100)

        Raises
        ------
        ValueError
            if camera=None is passed explicitly, instead pass an open 
            `firestudio.utils.camera_utils.Camera instance`

        Example usage
        -------------
```python
studio.set_ImageParams(
    this_setup_id='my_custom_setup',
    scale_bar=False,
    figure_label='high redshift')
```"""
        

        default_kwargs = {
            'camera':None,
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

        Returns
        ------- 
        abg_python.galaxy.metadata_utils.Metadata
            cache file for storing image maps

        Raises
        ------
        IOError
            if self.snapnum and self.sim_name are not set to disambiguate the cache file
        """
        

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

                self.metadata = self.cache_file = Metadata(
                    os.path.join(self.datadir,h5name),
                    loud=False) 
        else:
            raise IOError("Need to set self.snapnum and self.sim_name to disambiguate cache_file")

        return self.metadata

    def print_ImageParams(self):
        """ Prints the current image parameters."""
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
        """ stringifies image parameters and combines them to 'uniquely' 
            hash this combination of input in the cache file

        Returns
        -------
        str
            a stringified combination of image parameters
        """

        camera:Camera = self.camera
        ## uniquely identify this projection setup using a simple "hash"
        camera_string = 'camera'
        for vector in [camera.camera_pos,camera.camera_focus,camera.camera_north]:
            for component in vector: camera_string += f"_{component:0.2f}"
        
        self.this_setup_id = (
        "npix%d_half_width%.2fkpc_zmin%.2fkpc_zmax%.2fkpc_%s_aspect%.2f"%(
            self.pixels, 
                np.round(camera.frame_half_width,decimals=2),
                np.round(camera.zmin,decimals=2),
                np.round(camera.zmax,decimals=2),
                camera_string,
                np.round(self.aspect_ratio,decimals=2)))

        return self.this_setup_id

    def computeFrameBoundaries(self):
        """ Uses the camera to calculate the minimum and maximum 
            x, y, and z limits as well as the physical resolution of the image.
            

            Attributes set
            --------------
            self.Xmin, self.Xmax -- 
            self.Ymin, self.Ymax -- 
            self.Zmin, self.Zmax -- 
            self.npix_x, self.npix_y -- 
            self.Acell -- 
                
        """
        ## +- camera_dist limits -> 45 degree FOV
        self.Xmin,self.Xmax = -self.camera.frame_half_width,self.camera.frame_half_width

        self.Ymin,self.Ymax = np.array(
            [-self.camera.frame_half_width,self.camera.frame_half_width])*self.aspect_ratio

        self.Zmin,self.Zmax = -self.camera.zmin,self.camera.zmax

        ## Set image size 
        self.npix_x   = self.pixels #1200 by default
        self.npix_y   = int(self.pixels*self.aspect_ratio) #1200 by default

        self.Acell = (self.Xmax-self.Xmin)/self.npix_x * (self.Ymax-self.Ymin)/self.npix_y