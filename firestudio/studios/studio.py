from __future__ import print_function
import os
import numpy as np 
import h5py

import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec

from abg_python.snapshot_utils import openSnapshot
from abg_python.cosmo_utils import load_AHF
from abg_python.all_utils import append_function_docstring

from abg_python.galaxy.gal_utils import Galaxy
from abg_python.galaxy.metadata_utils import Metadata,metadata_cache

from firestudio.utils.stellar_utils.load_stellar_hsml import get_particle_hsml

class Drawer(object):
    
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

        ## will check relevant flags internally
        self.addText(ax)

        ## turn off the axis if asked
        if self.noaxis:
            ax.axis('off')

####### image utilities #######
    def addScaleBar(self,image):

        image_length =2*self.frame_half_width # kpc
        ## set scale bar length
        if image_length > 15 : 
            scale_line_length = 5
            self.scale_label_text = r"$\mathbf{5 \, \rm{kpc}}$"

        elif image_length > 1.5 : 
            scale_line_length = 1.0 
            self.scale_label_text = r"$\mathbf{1 \, \rm{kpc}}$"

        else:
            scale_line_length = .1
            self.scale_label_text = r"$\mathbf{100 \, \rm{pc}}$"

        # Convert to pixel space
        length_per_pixel = (self.Xmax - self.Xmin) / self.npix_x
        scale_line_length_px = int(scale_line_length / length_per_pixel)

        # Position in terms of image array indices
        scale_line_x_start = int(0.05 * self.npix_x)
        scale_line_x_end = min(scale_line_x_start + scale_line_length_px,self.npix_x)
        scale_line_y = int(0.02 * self.npix_y)

        # Go through pixels for scale bar, setting them to white
        for x_index in range(scale_line_x_start, scale_line_x_end):
            image[scale_line_y:scale_line_y+6, x_index,:3] = 1
        return image

    def addText(self,ax):
        ## handle any text additions
        if self.figure_label is not None:
        ## plot the  figure label in the top right corner
            label = ax.text(
                0.95, 0.92,
                self.figure_label,
                fontsize = self.fontsize,
                transform = ax.transAxes,
                ha='right')
            label.set_color('white')

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
            label2.set_color('white')

    def renormalizeTransposeImage(self,image,min_val,max_val,quantity_name):
        print('min_%s = '%quantity_name,min_val)
        print('max_%s = '%quantity_name,max_val)

        print('Image range (%s): '%quantity_name,np.min(image),np.max(image))
        image = image - min_val
        image = image / (max_val - min_val)
        
        ## clip anything outside the range
        image[image < 0.0] = 0.0
        image[image > 1.0] = 1.0
        image = image*255.0

        print('Image range (8bit): ',np.min(image),np.max(image))

        ## cast to integer to use as indices for cmap array
        image = image.astype(np.uint16) 

        return image.T

    def __saveFigure(
        self,
        ax,
        image_name,
        ):

        ## save the figure if asked
        savefig_args={} 

        if self.noaxis:
            ## remove whitespace around the axis, apparently the x/y origin is offset in pixel 
            ## space and so the bounding box doesn't actually reflect the left/bottom edge of the 
            ## axis
            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.yaxis.set_major_locator(plt.NullLocator())
            savefig_args['bbox_inches']='tight'
            savefig_args['pad_inches']=0

        image_name = "%s_%03d_%dkpc.png" % (image_name,self.snapnum, 2*self.frame_half_width)

        ax.get_figure().savefig(
            os.path.join(self.image_dir,image_name),dpi=300,
            **savefig_args)

class Studio(Drawer):
    """`FIREstudio` parent class that regularizes image setup, rotation, 
        caching, etc between `GasStudio` and `StarStudio` classes. 
    """

    def __init__(
        self,
        datadir, ## directory to put intermediate and output files 
        snapnum, ## snapshot number 
        cache_file_name = None, ##  the name of the file to save maps to
        gas_snapdict=None, ## open galaxy gas snapshot dictionary 
        star_snapdict=None, ## open galaxy star snapshot dictionary
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

star_snapdict['Coordinates'] ## coordinates of the particles
star_snapdict['Metallicity'] ## metallicity (mass fractions) of the particles 
star_snapdict['Masses']
star_snapdict['AgeGyr'] ## masses of the particles in 1e10 solar masses
```

(and ideally `'SmoothingLengths'`, in the same units as coordinates for both, but these can be calculated). 

            Input:

                datadir - directory to put intermediate and output files 
                snapnum - snapshot number (feel free to lie if you aren't using
                    FIRE_studio to open a snapshot, it is needed for cache file name though)
                cache_file_name=None, ##  the name of the file to save maps to
                    defaults to proj_maps_%03d.hdf5, changes when you use set_DataParams
                    to choose a snapshot.
                gas_snapdict=None - an open galaxy gas snapshot dictionary 
                star_snapdict=None - an open galaxy star snapshot dictionary 

                **set_ImageParams_kwargs -- keyword arguments that will be passed to
                    set_ImageParams

            Output:

                None"""
        
        if cache_file_name is None:
            cache_file_name = "proj_maps"

        self.cache_file_name = cache_file_name
        self.cache_file = None
        
        ## IO stuff
        self.datadir = datadir

        ## make the datadirectory if it doesn't exist
        if not os.path.isdir(self.datadir):
            os.makedirs(self.datadir)

        self.gas_snapdict = gas_snapdict
        self.star_snapdict = star_snapdict

        ## create, if necessary, directories to store intermediate and output files,
        ##  this could get crowded! sets self.image_dir and self.projection_dir
        #self.makeOutputDirectories(datadir)

        ## initialize the object with some default image params that will
        ##  make a face-on image of a galaxy, can always set these manually
        ##  with  external calls to set_ImageParams
        self.set_ImageParams(use_defaults=True,snapnum=snapnum,**kwargs)

####### I/O functions #######
    def load_SnapshotData(
        self,
        gas_mask=None,
        star_mask=None,
        **kwargs):
        """Binds simulation output to self.gas_snapdict and self.star_snapdict.

            Input:

                None

            Output:

                None"""

        ## determine if we need to open any snapshot data
        if (self.gas_snapdict is None or 
            self.gas_snapdict['snapnum'] != snapnum ): ## haven't loaded this data yet, or we are replacing it
            self.__get_snapdicts(
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
            self.masked_star_snapdict = filterDictionary(self.star_snapdict,mask)

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
            ['%s_SmoothingLengths'],
            use_metadata=use_metadata,
            save_meta=save_meta,
            assert_cached=assert_cached,
            loud=loud)
        def compute_HSML(self):
            snapdict = getattr(self,snapdict_name+'_snapdict')
            pos = snapdict['Coordinates']
            smoothing_lengths = get_particle_hsml(pos[:,0],pos[:,1],pos[:,2])
            return smoothing_lengths

        return_value = compute_HSML(self,**kwargs)
        ## TODO why am i getting back tuples...?
        if type(return_value == tuple):
            return_value  = return_value[0]
        return return_value

    def __get_snapdicts(
        self,
        sim_name,
        snapdir,snapnum,
        load_stars = 0,
        keys_to_extract = None,
        star_keys_to_extract = None,
        use_saved_subsnapshots=False,
        del_stars = False,
        **kwargs):

        galaxy = Galaxy(
            sim_name,
            snapdir,
            snapnum,
            datadir=self.datadir, ## not the same as where the FIREstudio cache is, I think... TODO
            loud_metadata=False, ## shh don't let them know
            save_header_to_table=False,
            **kwargs)

        ## handles opening the snapshot, centering it, and rotating it to be face-on.
        galaxy.extractMainHalo(
            use_saved_subsnapshots=use_saved_subsnapshots,
            save_meta=False,
            **kwargs) ## metadata cache will pull only the good keys out

        ## bind the snapshot dictionaries
        self.gas_snapdict = galaxy.sub_snap
        self.star_snapdict = galaxy.sub_star_snap

        ## just get rid of it now that we've opened it
        del galaxy

    def set_ImageParams(
        self,
        this_setup_id=None,
        use_defaults=False,
        edgeon=False,
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
                edgeon = False -- currently not implemented
                loud = True -- 

                frame_half_width = 15 --  half-width of image in x direction
                frame_depth = 15 -- z-depth of image (thickness is 2*frame_depth)
                frame_center = None -- center of frame in data space

                theta = 0,phi = 0,psi = 0 -- euler rotation angles in degrees

                aspect_ratio = 1 -- shape of image, y/x TODO figure out if this is necessary to pass?
                pixels = 1200 -- pixels in x direction, resolution of image
                figure_label = '' -- string to be put in upper right corner
                scale_bar = True -- flag to plot length scale bar in lower left corner
                noaxis = True -- turns off axis ticks
                savefig = None -- save the image as a png if passed a string
                fontsize = 12 
                
                snapdir - path to simulation output
                snapnum -  which snapshot to open
                sim_name - name of simulation (i.e. m12i_res7100)

            Output:

                None

Example usage:
```python
studio.set_ImageParams(
    this_setup_id='my_custom_setup',
    theta=90,
    scale_bar=False,
    figure_label='high redshift')
```"""

        default_kwargs = {
            'frame_half_width':15, ## half-width of image in x direction
            'frame_depth':None, ## z-depth of image (thickness is 2*frame_depth)
            'frame_center':None, ## center of frame in data space
            'theta':0,'phi':0,'psi':0, ## euler rotation angles
            'aspect_ratio':1, ## shape of image, y/x TODO figure out if this is necessary to pass?
            'pixels':1200, ## pixels in x direction, resolution of image
            'figure_label':'', ## string to be put in upper right corner
            'scale_bar':True,  ## flag to plot length scale bar in lower left corner
            'noaxis':True, ## turns off axis ticks
            'savefig':None, ## save the image as a png if passed a string
            'fontsize':12,  ## font size of scale bar and figure label
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
                if loud:
                    print("setting",kwarg,
                        'to user value of:',value)
                ## set it to the object
                setattr(self,kwarg,value)

        if use_defaults:
            ## if we haven't already set frame center by passed kwarg
            ##  need to replace None with [0,0,0]
            if ('frame_center' in default_kwargs and 
                default_kwargs['frame_center'] is None):
                default_kwargs['frame_center'] = np.zeros(3)

            ##  need to replace frame_depth with frame_half_width
            if ('frame_depth' in default_kwargs and 
                default_kwargs['frame_depth'] is None):

                ## take the default value
                if 'frame_center' in default_kwargs:
                    default_kwargs['frame_depth'] = default_kwargs['frame_half_width']
                ## take the value that was passed and set above
                else:
                    default_kwargs['frame_depth'] = self.frame_half_width
             
            ## set the remaining image parameters to their default values
            for default_arg in default_kwargs:
                value = default_kwargs[default_arg]
                if loud:
                    print("setting",default_arg,
                        'to default value of:',value)
                setattr(self,default_arg,value)


        if edgeon:
            raise NotImplementedError("need to fix this")
            print("Drawing an edgeon view, rotating theta = 90")
            self.theta+=90
            self.aspect_ratio*=0.5
            self.__identifyThisSetup()
            self.computeFrameBoundaries()
         
        ## determine the edges of our frame so we can cull the rest later,
        ##  also set cell area info
        self.computeFrameBoundaries()
        
        ## identify the combination of these parameters
        if this_setup_id is None:
            this_setup_id = self.__identifyThisSetup()

        self.this_setup_id = this_setup_id

        self.set_CacheFile()

    def set_CacheFile(self):
        """ Creates the cache hdf5 file. Requires self.snapnum and sim_name be set.

            Input:

                None

            Output:

                cache_file -- abg_python.galaxy.metadata_utils.Metadata object"""

        if self.snapnum is not None:
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
            raise IOError("Need to set self.snapnum to disambiguate cache_file")

        return self.metadata

    def print_ImageParams():
        """ Prints current image setup to console.

            Input:

                None

            Output:

                None"""
        default_kwargs = [
            'frame_half_width'
            'frame_depth'
            'frame_center'
            'theta','phi','psi',
            'aspect_ratio',
            'pixels', 
            'figure_label', 
            'scale_bar',  
            'noaxis'
            'savefig'
            'fontsize'
            'snapdir',
            'snapnum']

        ## print the current value, not the default value
        for arg in default_kwargs:
            print(arg,'=',getattr(self,arg))

    def __identifyThisSetup(self):
        ## uniquely identify this projection setup using a simple "hash"
        self.this_setup_id = (
	"npix%d_width%.2fkpc_depth%.2fkpc_x%.2f_y%.2f_z%.2f_theta%.2f_phi%.2f_psi%.2f_aspect%.2f"%(
	    self.pixels, 
            np.round(2*self.frame_half_width,decimals=2),
            np.round(self.frame_depth,decimals=2),
	    np.round(self.frame_center[0],decimals=2),
            np.round(self.frame_center[1],decimals=2),
            np.round(self.frame_center[2],decimals=2),
	    np.round(self.theta,decimals=2),
            np.round(self.phi,decimals=2),
            np.round(self.psi,decimals=2),
            np.round(self.aspect_ratio,decimals=2)
            ))
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
        self.Xmin,self.Xmax = self.frame_center[0] + np.array(
            [-self.frame_half_width,self.frame_half_width])

        self.Ymin,self.Ymax = self.frame_center[1] + np.array(
            [-self.frame_half_width,self.frame_half_width])*self.aspect_ratio

        self.Zmin,self.Zmax = self.frame_center[2] + np.array(
            [-self.frame_depth,self.frame_depth])

        ## Set image size 
        self.npix_x   = self.pixels #1200 by default
        self.npix_y   = int(self.pixels*self.aspect_ratio) #1200 by default

        self.Acell = (self.Xmax-self.Xmin)/self.npix_x * (self.Ymax-self.Ymin)/self.npix_y

    def cullFrameIndices(
        self,
        Coordinates):

        print("TODO:Need to decide if  we want to rotate before or after culling...")
        ## extract a cube of particles that are in relevant area
        ind_box = ((Coordinates[:,0] > self.Xmin) & (Coordinates[:,0] < self.Xmax) &
                   (Coordinates[:,1] > self.Ymin) & (Coordinates[:,1] < self.Ymax) &
                   (Coordinates[:,2] > self.Zmin) & (Coordinates[:,2] < self.Zmax))

        return ind_box

    def rotateEuler(self,theta,phi,psi,pos):
        pos-=self.frame_center
        ## if need to rotate at all really -__-
        if theta==0 and phi==0 and psi==0:
            return pos
        # rotate particles by angle derived from frame number
        pi        = 3.14159265
        theta_rad = pi*theta/ 1.8e2
        phi_rad   = pi*phi  / 1.8e2
        psi_rad   = pi*psi  / 1.8e2

        # construct rotation matrix
        rot_matrix = np.array([
            [np.cos(phi_rad)*np.cos(psi_rad), #xx
                -np.cos(phi_rad)*np.sin(psi_rad), #xy
                np.sin(phi_rad)], #xz
            [np.cos(theta_rad)*np.sin(psi_rad) + np.sin(theta_rad)*np.sin(phi_rad)*np.cos(psi_rad),#yx
                np.cos(theta_rad)*np.cos(psi_rad) - np.sin(theta_rad)*np.sin(phi_rad)*np.sin(psi_rad),#yy
                -np.sin(theta_rad)*np.cos(phi_rad)],#yz
            [np.sin(theta_rad)*np.sin(psi_rad) - np.cos(theta_rad)*np.sin(phi_rad)*np.cos(psi_rad),#zx
                np.sin(theta_rad)*np.cos(psi_rad) - np.cos(theta_rad)*np.sin(phi_rad)*np.sin(psi_rad),#zy
                np.cos(theta_rad)*np.cos(phi_rad)]#zz
            ]).astype(np.float32)

        n_box = pos.shape[0]

        ## rotate about each axis with a matrix operation
        pos_rot = np.matmul(rot_matrix,pos.T).T

        ## on 11/23/2018 (the day after thanksgiving) I discovered that 
        ##  numpy will change to column major order or something if you
        ##  take the transpose of a transpose, as above. Try commenting out
        ##  this line and see what garbage you get. ridiculous.
        pos_rot = np.array(pos_rot,order='C')
        
        ## add the frame_center back
        pos_rot+=self.frame_center

        ## can never be too careful that we're float32
        return pos_rot.astype(np.float32)

## append method docstrings to class docstring
append_function_docstring(Studio,Studio.__init__)
append_function_docstring(Studio,Studio.load_SnapshotData)
append_function_docstring(Studio,Studio.print_ImageParams)
append_function_docstring(Studio,Studio.computeFrameBoundaries)
append_function_docstring(Studio,Studio.get_HSML)
