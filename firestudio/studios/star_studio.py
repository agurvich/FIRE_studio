## builtin imports
import os
import sys 
import h5py
import matplotlib 
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import numpy as np 
import ctypes

## abg_python imports
from abg_python.plot_utils import addColorbar

## firestudio imports
from firestudio.studios.studio import Studio

from firestudio.utils.stellar_utils import raytrace_projection,load_stellar_hsml
import firestudio.utils.stellar_utils.make_threeband_image as makethreepic


class StarStudio(Studio):
    """
    Input:
        snapdir - location that the snapshots live in
        snapnum - snapshot number

        frame_half_width - half-width of image in data space
        frame_depth - half-depth of image in data space 

    Optional:
        overwrite=0 - flag to overwrite the intermediate grid files
        datadir=None - directory to output the the intermediate grid files and output png

        maxden=1e-2 - controls the saturation of the image in a non-obvious way
        dynrange=100.0 - controls the saturation of the image in a non-obvious way
        color_scheme_nasa=True - flag to use nasa colors (vs. SDSS if false) 

        snapdict=None - Dictionary-like holding gas snapshot data, open from disk if None
        star_snapdict=None - Dictionary-like holding star snapshot data, open from disk if None
        intermediate_file_name="proj_maps" ##  the name of the file to save maps to
    """ + "------- Studio\n" + Studio.__doc__

    def __init__(
        self,
        snapdir,snapnum, ## snapshot directory and snapshot number
        datadir, ## directory to put intermediate and output files
        frame_half_width, ## half-width of image in x direction
        frame_depth, ## z-depth of image (thickness is 2*frame_depth)
        maxden = 1.0e-2, ## controls the saturation of the image in a non-obvious way
        dynrange = 100.0, ## controls the saturation of the image in a non-obvious way
        color_scheme_nasa = True, ## flag to use nasa colors (vs. SDSS if false)
        star_snapdict = None, ## provide an open snapshot dictionary to save time opening
        snapdict = None, ## provide an open snapshot dictionary to save time opening
        **kwargs):

        self.maxden = maxden
        self.dynrange = dynrange
        self.color_scheme_nasa = color_scheme_nasa

        self.star_snapdict=star_snapdict
        self.snapdict=snapdict

        ## call Studio's init
        super().__init__(
            snapdir,snapnum,
            datadir,
            frame_half_width,
            frame_depth,
            **kwargs)

####### makeOutputDirectories implementation #######
    def makeOutputDirectories(self,datadir):
        print('Drawing %s:%d'%(self.snapdir,self.snapnum)+' to:%s'%self.datadir)

        ## make the path to store all plots in the datadir
        path = os.path.join(datadir,'Plots')
        if not os.path.isdir(path):
            os.mkdir(path)

        ## make the path to store final images
        self.image_dir = os.path.join(path,'HubbleImages')
        if not os.path.isdir(self.image_dir):
            os.mkdir(self.image_dir)

        ## make the paths to store projection hdf5 files
        self.projection_dir = os.path.join(path,'Projections')
        if not os.path.isdir(self.projection_dir):
            os.mkdir(self.projection_dir)

####### projectImage implementation #######
    def projectImage(self,image_names):

        ## open snapshot data if necessary
        if self.snapdict is None or self.star_snapdict is None:
            ## open and bind them in Studio's openSnapshot
            self.openSnapshot(
                load_stars = True,
                keys_to_extract = 
                    ['Coordinates',
                    'Masses',
                    'SmoothingLength', 
                    'Metallicity']+
                    ['Velocities','Density']*self.extract_galaxy,
                star_keys_to_extract = 
                    ['Coordinates',
                    'Masses',
                    'AgeGyr',
                    'Metallicity']+
                    ['Velocities']*self.extract_galaxy,
                )

        ## cull the particles outside the frame and cast to float32
        star_ind_box = self.cullFrameIndices(self.star_snapdict['Coordinates'])
        print(np.sum(star_ind_box),'many stars in volume')

        ## unpack the star information
        ## dont' filter star positions just yet
        star_pos = self.star_snapdict['Coordinates']

        ## try opening the stellar smoothing lengths, if we fail
        ##  let's calculate them and save them to the projection 
        ##  file
        try:
            with h5py.File(self.projection_file, "r") as handle:
                group = handle['PartType4']
                h_star = group['h_star'].value

            ## attempt to pass these indices along
            h_star = h_star[star_ind_box]
        except (KeyError,OSError,IndexError):
            print("Haven't computed stellar smoothing lengths...")
            h_star = load_stellar_hsml.get_particle_hsml(
                star_pos[:,0],star_pos[:,1],star_pos[:,2])

            ## write the output to an .hdf5 file
            with h5py.File(self.projection_file, "a") as handle:
                ## find a nice home in the hdf5 file for it
                if 'PartType4' not in handle.keys():
                    group = handle.create_group('PartType4')
                else:
                    group = handle['PartType4']

                ## overwrite existing h_star
                if 'h_star' in group.keys():
                    del group['h_star']
                group['h_star'] = h_star 

            h_star = h_star[star_ind_box]
            print("Done!")

        ## and now filter the positions
        star_pos = star_pos[star_ind_box]

        ## rotate by euler angles if necessary
        star_pos = self.rotateEuler(self.theta,self.phi,self.psi,star_pos)

        mstar =self.star_snapdict['Masses'][star_ind_box]
        ages = self.star_snapdict['AgeGyr'][star_ind_box]
        metals = self.star_snapdict['Metallicity'][:,0][star_ind_box]

        ## cull the particles outside the frame and cast to float32
        gas_ind_box = self.cullFrameIndices(self.snapdict['Coordinates'])
        print(np.sum(gas_ind_box),'many gas in volume')

        ## unpack the gas information
        gas_pos = self.snapdict['Coordinates'][gas_ind_box]

        ## rotate by euler angles if necessary
        gas_pos = self.rotateEuler(self.theta,self.phi,self.psi,gas_pos)

        mgas = self.snapdict['Masses'][gas_ind_box]
        gas_metals = self.snapdict['Metallicity'][:,0][gas_ind_box]
        h_gas = self.snapdict['SmoothingLength'][gas_ind_box]

        ## do the actual raytracing
        gas_out,out_u,out_g,out_r = raytrace_ugr_attenuation(
            star_pos[:,0],star_pos[:,1],star_pos[:,2],
            mstar,ages,metals,
            h_star,
            gas_pos[:,0],gas_pos[:,1],gas_pos[:,2],
            mgas,gas_metals,h_gas,
            pixels=self.pixels)

        ## write the output to an .hdf5 file
        self.writeImageGrid(
            out_u,'out_u',
            overwrite=self.overwrite)

        self.writeImageGrid(
            out_g,'out_g',
            overwrite=self.overwrite)

        self.writeImageGrid(
            out_r,'out_r',
            overwrite=self.overwrite)

####### produceImage implementation #######
    def produceImage(self,image_names):
        with h5py.File(self.projection_file, "r") as handle:
            this_group=handle[self.this_setup_id]
            out_u = this_group['out_u'].value
            out_g = this_group['out_g'].value
            out_r = this_group['out_r'].value

        ## open the hdf5 file and load the maps
        image24, massmap = makethreepic.make_threeband_image_process_bandmaps(
            out_r,out_g,out_u,
            maxden=self.maxden,
            dynrange=self.dynrange,
            pixels=self.pixels,
            color_scheme_nasa=self.color_scheme_nasa,
            color_scheme_sdss=not self.color_scheme_nasa)

        ## for some reason it's rotated 90 degrees...? kind of like transposed but different
        ##  need to take the rows of the output image and make them the columns, iteratively,
        ##  for now... 
        #image24=np.rot90(image24,k=1,axes=(0,1))
        final_image = np.transpose(image24,axes=(1,0,2))
            
        self.final_image = final_image

        return final_image

##### Image projection stuff
## Stellar light attenuation projection
def raytrace_ugr_attenuation(
    x,y,z,
    mstar,ages,metals,
    h_star, 
    gx,gy,gz,
    mgas, gas_metals,
    h_gas,
    pixels = 1200,
    xlim = None, ylim = None, zlim = None
    ):

    ## setup boundaries to cut-out gas particles that lay outside
    ## range
    if xlim is None:
        xlim = [np.min(x),np.max(x)]
    if ylim is None:
        ylim = [np.min(y),np.max(y)]
    if zlim is None:
        zlim = [np.min(z),np.max(z)]

#   band=BAND_ID; # default=bolometric
#   j = [  0,  6,  7,  8,  9, 10, 11, 12, 13,  1,   2,   3,   4,   5] # ordering I'm used to
#   i = [  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10,  11,  12,  13] # ordering of this
#   band_standardordering = band
#   band = j[band]
#   if (band > 13): 
#       print 'BAND_ID must be < 13'; 
#       return 0;
#   
#   b=['Bolometric', \
#   'Sloan u','Sloan g','Sloan r','Sloan i','Sloan z', \
#   'Johnsons U','Johnsons B', 'Johnsons V','Johnsons R','Johnsons I', \
#   'Cousins J','Cousins H','Cousins K']
    ## pick color bands by their IDs, see above
    BAND_IDS=[9,10,11]
    #gas_out,out_u,out_g,out_r = stellar_raytrace(
    return raytrace_projection.stellar_raytrace(
        BAND_IDS,
        x,y,z,
        mstar,ages,metals,
        h_star,
        gx,gy,gz,
        mgas, gas_metals,
        h_gas,
        pixels=pixels,
        xlim=xlim,ylim=ylim,zlim=zlim,
        ADD_BASE_METALLICITY=0.001*0.02,ADD_BASE_AGE=0.0003,
        IMF_SALPETER=0,IMF_CHABRIER=1
    )
