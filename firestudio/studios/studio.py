from __future__ import print_function
import os
import numpy as np 
import h5py

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from abg_python.snapshot_utils import openSnapshot
from abg_python.cosmo_utils import load_AHF
from abg_python.cosmoExtractor import diskFilterDictionary

shared_kwargs = [
    'snapdir=', #--snapdir: place where snapshots live
    'snapstart=', #--snapstart : which snapshot to start the loop at
    'snapmax=', #--snapmax : which snapshot to end the loop at

    ## snapshot opening and extraction
    'extract_galaxy=', #--extract_galaxy=False : flag to use abg_python.cosmoExtractor to extract main halo
    'ahf_path=', #--ahf_path : path relative to snapdir where the halo files are stored

    ## intermediate projection file options
    'datadir=', #--datadir: place to output frames to
    'overwrite=', #--overwrite: flag to  overwrite the cached projection if it exists
    'this_setup_id=', ## None - string that defines the projection setup, None by default means
    'h5prefix=', ## '' - a string that you can prepend to the projection filename if desired
    'intermediate_file_name=', ## None, the name of the file to save maps to

    ## image orientation and properties
    'frame_half_width=', #--frame_half_width : half width of frame in kpc
    'frame_depth=', #--frame_depth : half depth of frame in kpc
    'theta=','phi=','psi=', #--theta,phi,psi : euler angles for rotation
    'edgeon=', #--edgeon : flag for sticking a 90 degree edge on rotation underneath 
    'aspect_ratio=', ## the 'shape' of the image (y/x)
    'pixels=', #--pixels : how many pixels in each direction to use, defaults to 1200

    ## image annotation
    'figure_label=', #--figure_label: text to put in the upper right hand corner
    'noaxis=', #--noaxis : flag for removing axis and whitespace for just the pretty part
    'fontsize=', ## None - fontsize for all text in frame
    'scale_bar=', ##1 - should you plot a scale bar in the bottom left corner

    ## parallel multiprocessing 
    'multiproc=', #--multiproc : how many processes should be run simultaneously, keep in mind memory constraints
]

class Studio(object):
    """ 
    Input:
        snapdir - location that the snapshots live in
        snapnum - snapshot number
        frame_center - origin of image in data space 
        frame_half_width - half-width of image in data space
        frame_depth - half-depth of image in data space 

    Optional:
        theta=0- euler rotation angle
        phi=0- euler rotation angle
        psi=0 - euler rotation angle
        aspect_ratio=1 - the 'shape' of the image (y/x)
        pixels=1200 - the resolution of image (pixels x pixels)
        h5prefix='' - a string that you can prepend to the projection filename if desired
        fontsize=None - fontsize for all text in frame
        figure_label=None - what string should you put in the top right corner? 
        scale_bar=1 - should you plot a scale bar in the bottom left corner
        overwrite=False - flag to overwrite intermediate maps in projection file
        this_setup_id=None - string that defines the projection setup, None by default means
            it defaults to a gross combination of frame params + angles
        noaxis=0 - flag to turn off axis (1=off 0=on)
        savefig=1 - flag to save figure to datadir (default snapdir, but can be a kwarg)
        ahf_path=None - path relative to snapdir where the halo files are stored
            defaults to snapdir/../halo/ahf
        extract_galaxy=False - flag to extract the main galaxy using abg_python.cosmoExtractor
        intermediate_file_name=None - the name of the file to save maps to
    """
    def __init__(
        self,
        snapdir,snapnum, ## snapshot directory and snapshot number
        datadir, ## directory to put intermediate and output files
        frame_half_width, ## half-width of image in x direction
        frame_depth, ## z-depth of image (thickness is 2*frame_depth)
        frame_center = None, ## center of frame in data space
        theta=0,phi=0,psi=0, ## euler rotation angles
        aspect_ratio = 1, ## shape of image, y/x
        pixels=1200, ## pixels in x direction
        h5prefix=None, ## string to prepend to projection file
        fontsize = 12,  ## font size of scale bar and figure label
        figure_label = None, ## string to be put in upper right corner
        scale_bar = True,  ## flag to plot length scale bar in lower left corner
        overwrite = False, ## flag to overwrite intermediate flags
        this_setup_id = None, ## string identifier in the intermediate projection file
        noaxis = True, ## turns off axis ticks
        savefig = True, ## save the image as a png
        ahf_path = None, ## path relative to snapdir where the halo files are stored
        extract_galaxy = False, ## uses halo center to extract region around main halo
        intermediate_file_name = None, ##  the name of the file to save maps to
        **kwargs
        ):
        
        ## handle defaults
        if frame_center is None:
            frame_center = np.zeros(3)

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
            self.addScaleBar(self.final_image)

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

        ## turn off the axis if asked
        if self.noaxis:
            ax.axis('off')

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

            else:
                ## appending another image, or overwriting a map, cool!
                this_group = h5file[self.this_setup_id]

                ## might want to overwrite a quantity map
                if image_name in this_group.keys():
                    if overwrite:
                        del this_group[image_name]
                    else:
                        if not np.all(image == this_group[image_name]):
                            raise IOError(
                            "%s already exists with overwrite=False."%
                            image_name)
                        else:
                            return


                ## save this new quantity
                this_group[image_name] = image 

    def saveFigure(
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

####### data utilities #######
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

#### FUNCTIONS THAT SHOULD BE OVERWRITTEN IN SUBCLASSES
    def makeOutputDirectories(self,datadir):
        raise NotImplementedError("Studio is a base-class and this method must be implemented in a child.")

    def projectImage(self):
        raise NotImplementedError("Studio is a base-class and this method must be implemented in a child.")

    def produceImage(self):
        raise NotImplementedError("Studio is a base-class and this method must be implemented in a child.")

