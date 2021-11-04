from __future__ import print_function
import os
import numpy as np 
import h5py

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from abg_python.snapshot_utils import openSnapshot
from abg_python.cosmo_utils import load_AHF
from abg_python import append_function_docstring,getThetasTaitBryan,filterDictionary
from abg_python.plot_utils import nameAxes

        if intermediate_file_name is None:
            intermediate_file_name = "proj_maps"
        
        if h5prefix is None:
            h5prefix = ''

        if figure_label is None:
            figure_label = ''

        print("extra kwargs:\n",list(kwargs.keys()))

        ## IO stuff
        self.snapdir = snapdir
        self.snapnum = snapnum
        self.datadir = datadir

        ## bind this setup's parameters
        self.frame_center = frame_center
        self.frame_half_width = frame_half_width
        self.frame_depth = frame_depth
        self.theta,self.phi,self.psi = theta,phi,psi
        self.aspect_ratio = aspect_ratio
        self.pixels = pixels

        ## plot/overlay parameters
        self.fontsize = fontsize
        self.figure_label = figure_label
        self.scale_bar = scale_bar
        self.scale_bar_length = scale_bar_length

        ## should we overwrite intermediate files?
        self.overwrite = overwrite

        self.noaxis = noaxis
        self.savefig = savefig
        self.ahf_path = ahf_path
        self.extract_galaxy = extract_galaxy

        ## create, if necessary, directories to store intermediate and output files,
        ##  this could get crowded! sets self.image_dir and self.projection_dir
        self.makeOutputDirectories(datadir)

        h5name=h5prefix+intermediate_file_name+"_%03d.hdf5"% snapnum
        self.projection_file = os.path.join(self.projection_dir,h5name)

        ## determine the edges of our frame so we can cull the rest later
        self.computeFrameBoundaries()

        ## create the unique identifier string for this setup
        if this_setup_id is None:
            self.identifyThisSetup()
        else:
            self.this_setup_id = this_setup_id
         
    def renderFaceAppendEdgeViews(self,image_names):
        fig = plt.figure()
        gs = gridspec.GridSpec(3,1)

        axs = [
            fig.add_subplot(gs[0:2,0]),
            fig.add_subplot(gs[2:3,0])
        ]

        gs.update(wspace=0,hspace=-0.1) ## NOTE why doesn't hspace = 0 work???

        ## do the face on view
        self.scale_bar = 0
        self.render(
            axs[0],
            image_names)

        ## do the edge on view
        self.figure_label = None
        self.scale_bar = 1

        ## set figure size extended for the edgeon view
        fig.set_size_inches(6,9)
        self.render(
            axs[1],
            image_names,
            edgeon=True)

        plt.close(fig)

    def render(
        self,
        ax,
        image_name,
        edgeon=0,
        assert_cached=False):
        if ax is None:
            fig,ax = plt.figure(),plt.gca()
        else:
            fig = ax.get_figure()

        if edgeon:
            print("Drawing an edgeon view, rotating theta = 90")
            self.theta+=90
            self.aspect_ratio*=0.5
            self.identifyThisSetup()
            self.computeFrameBoundaries()
        
        ## unpack image_name if there are  multiple images that need
        ##  to be combined. self.produce_image should handle this invisibly
        if type(image_name) == list:
            ## the last string in the list will be the one used to name 
            ##  the actual .png at the end. so you could even do 
            ##  ['columDensityMap','frame'] if you wanted a simple
            ##  "1 color" image but wanted to change the name of the
            ##  .png
            image_names,image_name = image_name[:-1],image_name[-1]

        ## check if we've already projected this setup and saved it to intermediate file
        this_setup_in_projection_file = self.checkProjectionFile(image_names)

        ## allow the user to require that the setup is cached
        if assert_cached and not this_setup_in_projection_file:
            raise AssertionError("User required that setup was cached -- assert_cached=True")

        ## project the image using a C routine
        if not this_setup_in_projection_file or self.overwrite:
            self.projectImage(image_names)

        ## remap the C output to RGB space
        self.final_image = self.produceImage(image_names)

        ## plot that RGB image and overlay scale bars/text
        self.plotImage(ax,image_names)

        ## save the image
        if self.savefig:
            self.saveFigure(ax,image_name)

        ## return these to their previous values, because why not?
        if edgeon:
            self.theta-=90
            self.aspect_ratio*=2
            self.identifyThisSetup()
            self.computeFrameBoundaries()

    def plotImage(
        self,
        ax,
        image_names,
        cbar_label=None,
        **kwargs): 

        ## fill the pixels of the the scale bar with white
        if self.scale_bar:
            self.addScaleBar(self.final_image,self.scale_bar_length)

        ## main imshow call
        imgplot = ax.imshow(
            self.final_image, 
            extent = (self.Xmin,self.Xmax,self.Ymin,self.Ymax),
            origin = 'lower')

        ## set the y units to equal the x units-- make sure nothing
        ##  is being skewed. UNRELATED to self.aspect_ratio
        ax.set_aspect(1)

        ## will check relevant flags internally
        self.addText(ax)

####### image utilities #######
    def addScaleBar(self,image):

        image_length =2*self.frame_half_width # kpc
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

####### I/O functions #######
    def openSnapshot(
        self,
        load_stars = 0,
        keys_to_extract = None,
        star_keys_to_extract = None):

        if (self.snapdict is not None and
            ('star_snapdict' in self.__dict__ and self.star_snapdict is not None)):
            return
        elif not self.extract_galaxy:
            ## isolated galaxy huh? good choice. 
            ##  don't worry, cosmological will be overwritten in openSnapshot if 
            ##  HubbleParam != 1, none of the coordinates will be offset
            ##  and the galaxy won't be rotated or extracted
            snapdict = openSnapshot(
                self.snapdir,self.snapnum,
                ptype=0,cosmological=0,
                keys_to_extract=keys_to_extract)
            if load_stars:
                star_snapdict = openSnapshot(
                    self.snapdir,self.snapnum,
                    ptype=4,cosmological=0,
                    keys_to_extract=star_keys_to_extract)

                ## could just be a sub-snapshot that's been pre-extracted
                if not star_snapdict['cosmological']:
                    raise Exception("Need to open type 1 and type 2 parts")
        else:
            snapdict = openSnapshot(
                self.snapdir,self.snapnum,
                ptype=0,cosmological=1,
                keys_to_extract=keys_to_extract)

            if load_stars:
                star_snapdict = openSnapshot(
                    self.snapdir,self.snapnum,
                    ptype=4,cosmological=1,
                    keys_to_extract=star_keys_to_extract)
            else:
                star_snapdict = None

            ## cosmological snapshot it is then... 
            scom,rvir,vesc = load_AHF(
                self.snapdir,self.snapnum,
                snapdict['Redshift'],
                ahf_path=self.ahf_path)

            ## filter all the keys in the snapdict as necessary to extract a spherical volume
            ##  centered on scom (the halo center), using 3*frame_half_width
            diskFilterDictionary(
                star_snapdict if load_stars else None, 
                snapdict,
                radius=3*self.frame_half_width, ## particles to mask
                orient_radius=3*self.frame_half_width, ## particles to orient on (in principle could orient on inner region and want larger region)
                scom=scom,
                orient_stars = load_stars)

        ## bind the snapdicts
        self.snapdict = snapdict
        if load_stars:
            self.star_snapdict = star_snapdict

    def identifyThisSetup(self):
        ## uniquely identify this projection setup
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

    def checkProjectionFile(self,image_names):
        try:
            with h5py.File(self.projection_file,'r') as handle:
                for group in handle.keys():
                    this_group = handle[group]
                    flag = True
                    for key in ['npix_x','frame_half_width','frame_depth',
                        'frame_center','theta','phi','psi','aspect_ratio']:
                        variable = getattr(self,key)
                        if key not in this_group.keys():
                            flag = False
                            break

                        ## read the value in the hdf5 file and compare to variable
                        if key not in ['npix_x']:
                            ## key is not an integer/string so we have to round it somehow
                            flag = flag and np.all(
                                np.round(this_group[key].value,decimals=2) == np.round(variable,decimals=2))
                        else:
                            ## key can be directly compared
                            flag = flag and this_group[key].value == variable
                    ## found the setup one we wanted, does it have all the images we want?
                    if flag:
                        for image_name in image_names:
                            flag = flag and (image_name in this_group.keys())
                        return flag
            return 0 
        except IOError:
            return 0

    def writeImageGrid(
        self,
        image,
        image_name,
        overwrite=0):

        ## what should we call this setup? need a unique identifier
        ## let the user give it a
        ##	custom name through kwargs later on TODO

        with h5py.File(self.projection_file, "a") as h5file:
            if self.this_setup_id not in list(h5file.keys()):
                this_group = h5file.create_group(self.this_setup_id)
                ## save the maps themselves
                this_group[image_name]=image

                ## save the meta data
                this_group['npix_x']=self.npix_x
                this_group['frame_center']=self.frame_center
                this_group['frame_half_width']=self.frame_half_width
                this_group['frame_depth']=self.frame_depth
                this_group['theta']=self.theta
                this_group['phi']=self.phi
                this_group['psi']=self.psi
                this_group['aspect_ratio']=self.aspect_ratio
                ## TODO should I put in metadata that allows you to recreate
                ##  frames without access to the relevant snapshots?
                ##  e.g. current time/redshift

    def saveFigure(
        self,
        ax,
        image_name=None,
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

        
        if image_name is None:
            image_name = "%03d_%dkpc.pdf" % (self.snapnum, 2*self.frame_half_width)

        if 'png' not in image_name and 'pdf' not in image_name:
            image_name+='.pdf'

        ax.get_figure().savefig(
            os.path.join(self.datadir,image_name),dpi=300,
            **savefig_args)

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

                frame_half_width = 15 --  half-width of image in x direction
                frame_half_thickness = 15 -- half-thickness of image in z direction
                frame_center = None -- center of frame in data space

                theta = 0,phi = 0,psi = 0 -- euler rotation angles in degrees

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
    theta=90,
    scale_bar=False,
    figure_label='high redshift')
```"""

        default_kwargs = {
            'frame_half_width':15, ## half-width of image in x direction
            'frame_half_thickness':None, ## half-thickness of image in z direction
            'frame_center':None, ## center of frame in data space
            'theta':0,'phi':0,'psi':0, ## euler rotation angles
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
            ## if we haven't already set frame center by passed kwarg
            ##  need to replace None with [0,0,0]
            if ('frame_center' in default_kwargs and 
                default_kwargs['frame_center'] is None):
                default_kwargs['frame_center'] = np.zeros(3)

            ##  need to replace frame_half_thickness with frame_half_width
            if ('frame_half_thickness' in default_kwargs and 
                default_kwargs['frame_half_thickness'] is None):

                ## take the frame_half_width
                if 'frame_half_width' in default_kwargs:
                    default_kwargs['frame_half_thickness'] = default_kwargs['frame_half_width']
                ## take the value that was passed and set above
                else:
                    default_kwargs['frame_half_thickness'] = self.frame_half_width
             
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
            'frame_half_width':15, ## half-width of image in x direction
            'frame_half_thickness':None, ## half-thickness of image in z direction
            'frame_center':None, ## center of frame in data space
            'theta':0,'phi':0,'psi':0, ## euler rotation angles
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
            'frame_half_width',
            'frame_half_thickness',
            'frame_center',
            'theta','phi','psi',
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
        self.this_setup_id = (
	"npix%d_width%.2fkpc_depth%.2fkpc_x%.2f_y%.2f_z%.2f_theta%.2f_phi%.2f_psi%.2f_aspect%.2f"%(
	    self.pixels, 
            np.round(2*self.frame_half_width,decimals=2),
            np.round(self.frame_half_thickness,decimals=2),
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
        self.Xmin,self.Xmax = self.frame_center[0] + np.array(
            [-self.frame_half_width,self.frame_half_width])

        self.Ymin,self.Ymax = self.frame_center[1] + np.array(
            [-self.frame_half_width,self.frame_half_width])*self.aspect_ratio

        self.Zmin,self.Zmax = self.frame_center[2] + np.array(
            [-self.frame_depth,self.frame_depth])

        ## Set image size 
        self.npix_x   = self.pixels #1200 by default
        self.npix_y   = int(self.pixels*self.aspect_ratio) #1200 by default

    def cullFrameIndices(
        self,
        Coordinates):

        ## extract a cube of particles that are in relevant area
        print('extracting cube')
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

####### image utilities #######
    def addScaleBar(self,image,scale_bar_length):

        image_length =2*self.frame_half_width # kpc
        ## set scale bar length
        if scale_bar_length is None:
            if image_length > 15 : 
                scale_line_length = 5
                self.scale_label_text = r"$\mathbf{5 \, \rm{kpc}}$"

            elif image_length > 1.5 : 
                scale_line_length = 1.0 
                self.scale_label_text = r"$\mathbf{1 \, \rm{kpc}}$"

            else:
                scale_line_length = .1
                self.scale_label_text = r"$\mathbf{100 \, \rm{pc}}$"
        else:
            scale_line_length = scale_bar_length
            self.scale_label_text = (
                r"$\mathbf{" + 
                "{:.3g}".format( scale_bar_length ) +
                r" \, \rm{kpc}}$"
            )

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

#### FUNCTIONS THAT SHOULD BE OVERWRITTEN IN SUBCLASSES
    def makeOutputDirectories(self,datadir):
        raise NotImplementedError("Studio is a base-class and this method must be implemented in a child.")

    def projectImage(self):
        raise NotImplementedError("Studio is a base-class and this method must be implemented in a child.")

    def produceImage(self):
        raise NotImplementedError("Studio is a base-class and this method must be implemented in a child.")

