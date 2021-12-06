## builtin imports
import matplotlib 
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import numpy as np 
import copy

## abg_python imports
from abg_python.plot_utils import addColorbar,nameAxes
from abg_python import append_function_docstring,append_string_docstring,findIntersection
from abg_python.galaxy.metadata_utils import metadata_cache

## firestudio imports
from .studio import Studio

from ..utils.stellar_utils import (
    opacity_per_solar_metallicity, ## calculate dust opacity for arbitrary wavelengths
    read_band_lums_from_tables, ## calculate luminosities in specific tabulated bands
    stellar_raytrace, ## project along the LoS and attenuate by obscuring dust column for 3 bands simultaneously
    make_threeband_image_process_bandmaps) ## take 3 bands, overlay them, and map to RGB to make final image


class StarStudio(Studio):
    """ `FIREstudio` class for making mock hubble images with 
        attenuation along the line of sight. 

* [`StarStudio.get_mockHubbleImage`](#starstudioget_mockhubbleimage) 
* [`StarStudio.render`](#starstudiorender) 
* [`Studio.__init__`](#studio__init__) 
* [`Studio.set_ImageParams`](#studioset_imageparams)"""

    def __repr__(self):
        return 'StarStudio instance'

####### makeOutputDirectories implementation #######
    def set_ImageParams(
        self,
        use_defaults=False,
        loud=True,
        **kwargs):
        """Changes the parameters of the image. If `use_defaults=True` then 
            default values of the parameters will be set. Leave `use_defaults=False`
            to adjust only the keywords passed. 

            Input: 

                maxden = None --  controls the saturation of the image,
                    sets the upper limit of the "colorbar," defaults to 
                    the 99 %'ile of the image surface brightness
                dynrange = None  --  controls the saturation of the image,
                    sets the lower limit of the "colorbar" with respect to maxden,
                    defaults to the dynamic range between maxden and the 10th %'ile
                    of the image surface brightness
                color_scheme_nasa = True -- flag for switching between Hubble vs. SDSS images
                loud = True -- flagwhether print statements should show up on console.
            
            Output: 

                None

Example usage:
```python
starStudio.set_ImageParams(
    maxden=0.1,
    dynrange=10,
    figure_label='Hubble')
```"""

        default_kwargs = {
            'maxden' : None, ## 
            'dynrange' : None, ## controls the saturation of the image in a non-obvious way
            'color_scheme_nasa' : True} ## flag to use nasa colors (vs. SDSS if false)

        for kwarg in list(kwargs.keys()):
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
        super().set_ImageParams(use_defaults=use_defaults,**kwargs)

    append_function_docstring(set_ImageParams,Studio.set_ImageParams,prepend_string='passes `kwargs` to:\n')

    def print_ImageParams(self):
        """ Prints current image setup to console.

            Input:

                None

            Output:

                None"""

        default_kwargs = {
            'maxden' : 1.0e-2, ## 
            'dynrange' : 100.0, ## controls the saturation of the image in a non-obvious way
            'color_scheme_nasa' : True} ## flag to use nasa colors (vs. SDSS if false)

        ## print the current value, not the default value
        for arg in default_kwargs:
            print(arg,'=',getattr(self,arg))

        ## call the super class' print image params
        super().print_ImageParams()

    def quick_get_mockHubbleImage(
        self,
        use_metadata=True,
        save_meta=True,
        assert_cached=False,
        loud=True,
        lums=None,
        nu_effs=None,
        BAND_IDS=None
        ):

        if BAND_IDS is None:
            BAND_IDS=[1,2,3] ## used if corresponding column of lums is all 0s
       
        # apply filters, rotations, unpack snapshot data, etc...
        (kappas, lums,
            star_pos, mstar, ages, metals, h_star,
            gas_pos , mgas , gas_metals ,  h_gas) = self.prepareCoordinates(lums,nu_effs,BAND_IDS)

        star_xs,star_ys,star_zs = star_pos.T
        gas_xs,gas_ys,gas_zs = gas_pos.T

        xedges = np.linspace(self.Xmin,self.Xmax,self.npix_x,endpoint=True)
        yedges = np.linspace(self.Ymin,self.Ymax,self.npix_y,endpoint=True)

        ## area of cell
        dA = xedges[1] - xedges[0]
        dA*=dA

        outs = np.zeros((3,xedges.size-1,yedges.size-1))

        metal_mass_map,xedges,yedges = np.histogram2d(
            gas_xs,gas_ys,
            bins=[xedges,yedges],
            weights=mgas*gas_metals)

        for i,lum in enumerate(lums):
            outs[i],xedges,yedges = np.histogram2d(
                star_xs,star_ys,
                bins=[xedges,yedges],
                weights=lum)
        
            outs[i]*=np.exp(-kappas[i]*metal_mass_map/dA*5)

        unit_factor = 1e10/self.Acell

        return metal_mass_map*unit_factor,outs[0]*unit_factor,outs[1]*unit_factor,outs[2]*unit_factor

    def get_mockHubbleImage(
        self,
        use_metadata=True,
        save_meta=True,
        assert_cached=False,
        loud=True,
        **kwargs, 
        ):
        """Projects starlight and approximates attenuatation along line of sight
            into SDSS u, g, and r bands. 

            Input:

                use_metadata = True -- flag to search cache for result
                save_meta = True -- flag to cache the result
                assert_cached = False -- flag to require a cache hit
                loud = True -- whether cache hits/misses should be announced
                    to the console.

                lums = None -- manual input for luminosities
                nu_effs = None -- manual input for effective frequency of input luminosity band

                BAND_IDS = None -- index for luminosity band to project.
                    lambda_eff=np.array([
                        ## Bolometric (?)
                                1.e-5,  
                        ## SDSS u       g       r      i      z
                                3551. , 4686. , 6165., 7481., 8931.,
                        ##      U       B       V      R
                                3600. , 4400. , 5556., 6940., 
                        ##      I      J       H       K
                                8700., 12150., 16540., 21790.])

                    BAND_IDS[i] is ignored if nu_effs[i] is not None

            Output:

                gas_out -- total mass along LOS in pixel, in msun/kpc^2
                out_u -- total attenuated luminosity along LOS in pixel
                    in u band, in Lsun/kpc^2
                out_g -- total attenuated luminosity along LOS in pixel
                    in g band, in Lsun/kpc^2
                out_r -- total attenuated luminosity along LOS in pixel
                    in r band, in Lsun/kpc^2"""

        @metadata_cache(
            self.this_setup_id,  ## hdf5 file group name
            ['starMassesMap',
                'attenUMap', ## TODO naming it this could be confusing if BAND_IDS is different...
                'attenGMap',
                'attenRMap'],
            use_metadata=use_metadata,
            save_meta=save_meta,
            assert_cached=assert_cached,
            loud=loud,
            force_from_file=True)  ## read from cache file, not attribute of object
        def compute_mockHubbleImage(self,lums=None,nu_effs=None,BAND_IDS=None):

            if BAND_IDS is None:
                BAND_IDS=[1,2,3] ## used if corresponding column of lums is all 0s
            
            # apply filters, rotations, unpack snapshot data, etc...
            (kappas, lums,
                star_pos, mstar, ages, metals, h_star,
                gas_pos , mgas , gas_metals ,  h_gas) = self.prepareCoordinates(lums,nu_effs,BAND_IDS)

            ## do the actual raytracing
            gas_out,out_u,out_g,out_r = raytrace_ugr_attenuation(
                star_pos[:,0],star_pos[:,1],star_pos[:,2],
                mstar,ages,metals,
                h_star,
                gas_pos[:,0],gas_pos[:,1],gas_pos[:,2],
                mgas,gas_metals,h_gas,
                kappas,lums,
                pixels=self.pixels,
                QUIET=not self.master_loud,
                xlim = (self.Xmin, self.Xmax),
                ylim = (self.Ymin, self.Ymax),
                zlim = (self.Zmin, self.Zmax)
                )

            ## unit factor, output is in Lsun/kpc^2
            unit_factor = 1e10/self.Acell
            return gas_out*unit_factor, out_u*unit_factor, out_g*unit_factor, out_r*unit_factor
        return compute_mockHubbleImage(self,**kwargs)


    def prepareCoordinates(self,
        lums=None,
        nu_effs=None,
        BAND_IDS=None):
        
        ## unpack the star information
        ## dont' filter star positions just yet
        star_pos = self.star_snapdict['Coordinates']

        ## rotate by euler angles if necessary
        star_pos = self.camera.rotate_array(star_pos,offset=True)

        ## cull the particles outside the frame and cast to float32
        star_ind_box = self.cullFrameIndices(star_pos)
        if self.master_loud:
            print(np.sum(star_ind_box),'many star particles in volume')
        
        ## try opening the stellar smoothing lengths, if we fail
        ##  let's calculate them and save them to the projection 
        ##  file

        if "SmoothingLength" not in self.star_snapdict:
            Hsml = self.get_HSML('star')
            if Hsml.size != star_pos.shape[0]:
                Hsml = self.get_HSML('star',use_metadata=False,save_meta=True)
        else:
            Hsml = self.star_snapdict['SmoothingLength'] ## kpc
        ## attempt to pass these indices along
        h_star = Hsml[star_ind_box].astype(np.float32)

        ## and now filter the positions
        star_pos = star_pos[star_ind_box].astype(np.float32)

        mstar = self.star_snapdict['Masses'][star_ind_box].astype(np.float32)
        ages = self.star_snapdict['AgeGyr'][star_ind_box].astype(np.float32)
        metals = self.star_snapdict['Metallicity']
        if len(np.shape(metals)) > 1: metals = metals[:,0]
        metals = metals[star_ind_box].astype(np.float32)

        ## apply frame mask to band luminosities
        if lums is not None:
            lums = lums[:,star_ind_box]

        ## will fill any columns of lums with appropriate BAND_ID 
        ##  if nu_eff for that column is not None
        lums,nu_effs = read_band_lums_from_tables(
            BAND_IDS, 
            mstar,ages,metals,
            ## flag to return luminosity in each band requested without projecting
            nu_effs=nu_effs,
            lums=lums,
            QUIET=not self.master_loud)

        ## calculate the kappa in this band using:
        ##  Thompson scattering + 
        ##  Pei (1992) + -- 304 < lambda[Angstroms] < 2e7
        ##  Morrison & McCammon (1983) -- 1.2 < lambda[Angstroms] < 413
        ## get opacities and luminosities at frequencies we need:

        ## important to set KAPPA_UNITS appropriately: code loads opacity (kappa) for 
        ##   the bands of interest in cgs (cm^2/g), must be converted to match units of input 
        ##   mass and size. the default it to assume gadget units (M=10^10 M_sun, l=kpc)
        KAPPA_UNITS=2.08854068444 ## cm^2/g -> kpc^2/mcode
        kappas = [KAPPA_UNITS*opacity_per_solar_metallicity(nu_eff) for nu_eff in nu_effs]


        ## rotate by euler angles if necessary
        gas_pos = self.camera.rotate_array(self.gas_snapdict['Coordinates'],offset=True)

        ## cull the particles outside the frame and cast to float32
        gas_ind_box = self.cullFrameIndices(gas_pos)
        if self.master_loud:
            print(np.sum(gas_ind_box),'many gas particles in volume')

        ## unpack the gas information
        gas_pos = gas_pos[gas_ind_box].astype(np.float32)


        mgas = self.gas_snapdict['Masses'][gas_ind_box].astype(np.float32)
        gas_metals = self.gas_snapdict['Metallicity']
        if len(np.shape(gas_metals)) > 1: gas_metals = gas_metals[:,0]
        gas_metals = gas_metals[gas_ind_box].astype(np.float32)

        ## set metallicity of hot gas to 0 so there is no dust extinction
        temperatures = self.gas_snapdict['Temperature'][gas_ind_box]
        gas_metals[temperatures>1e5] = 0

        if "SmoothingLength" not in self.gas_snapdict:
            h_gas = self.get_HSML('gas')[gas_ind_box]
        else:
            h_gas = self.gas_snapdict['SmoothingLength'][gas_ind_box].astype(np.float32)

        return (kappas, lums,
                star_pos, mstar, ages, metals, h_star,
                gas_pos , mgas , gas_metals ,  h_gas)

####### produceImage implementation #######
    def render(
        self,
        ax=None,
        **kwargs):
        """Plots a mock hubble image, along with any annotations/scale bars,
            using the stored image parameters.

            Input: 

                ax = None -- axis to plot image to, if None will create a new figure
                quick = False -- flag to use a simple 2d histogram (for comparison or
                    for quick iteration as the user defines the image parameters)

            Output:

                ax -- the axis the image was plotted to
                final_image -- 2x2x3 RGB pixel array

Example usage:
```python
starStudio.render(plt.gca())
```"""

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
            self.saveFigure(ax,self.savefig)

        return ax,final_image

    def produceImage(
        self,
        quick=False,
        **kwargs):

        if not quick:
            gas_out,out_u,out_g,out_r = self.get_mockHubbleImage(**kwargs)
        else:
            gas_out,out_u,out_g,out_r = self.quick_get_mockHubbleImage(**kwargs)

        all_bands = np.concatenate([out_u,out_g,out_r])
        maxden_guess,dynrange_guess = self.predictParameters(all_bands=all_bands)
        if self.maxden is None:
            print("setting maxden to best guess %.2g"%maxden_guess)
            self.set_ImageParams(
                maxden=maxden_guess)

        if self.dynrange is None:
            print("setting dynrange to best guess %.2g"%dynrange_guess)
            self.set_ImageParams(
                dynrange=dynrange_guess)

        ## open the hdf5 file and load the maps
        image24, massmap = make_threeband_image_process_bandmaps(
            copy.copy(out_r),copy.copy(out_g),copy.copy(out_u),
            maxden=self.maxden,
            dynrange=self.dynrange,
            color_scheme_nasa=self.color_scheme_nasa,
            color_scheme_sdss=not self.color_scheme_nasa,
            QUIET=not self.master_loud)

        ## for some reason it's rotated 90 degrees...? kind of like transposed but different
        ##  need to take the rows of the output image and make them the columns, iteratively,
        ##  for now... 
        #image24=np.rot90(image24,k=1,axes=(0,1))
        final_image = np.transpose(image24,axes=(1,0,2))
        self.final_image = final_image

        return final_image

    def predictParameters(
        self,
        left_percentile=0.1,
        right_percentile=0.99,
        all_bands=None,
        ax=None,
        quick=False):
        """ Guesses what the "best" values for maxden and dynrange are from
            the distribution of surface brightnesses in the current image. 
            Looks for the left_percentile and right_percentile and returns
            right_percentile and the distance between it and left_percentile
            (in log space). 
            
            Input:

                left_percentile = 0.1 -- lower bound on image surface brightness percentile
                right_percentile = 0.99 --  upper bound on image surface brightness percentile
                ax = None -- optionally plots distribution of surface brightnesses
                    (in some units...) with overlay of percentiles and such.
                quick = False -- flag to use a simple 2d histogram (for comparison or
                    for quick iteration as the user defines the image parameters)
            
            Output:
                
                maxden -- maximum surface brightness of the image
                dynrange -- distance between maximum and minimum surface brightness
                     in log space. 
"""

        if (all_bands is None):
            ## read the luminosity maps
            if not quick:
                gas_out,out_u,out_g,out_r = self.get_mockHubbleImage()
            else:
                gas_out,out_u,out_g,out_r = self.quick_get_mockHubbleImage()

            all_bands = np.concatenate([out_u,out_g,out_r])

        ## concatenate the luminosity maps and take the log of the non-empty ones
        rats = np.log10(all_bands.flatten())
        rats = rats[np.isfinite(rats)]
        h,edges = np.histogram(rats,bins=1000)
    
        ## take the CDF to find left and right percentiles
        cumulative = np.array(np.cumsum(h))
        cumulative=cumulative/cumulative[-1]

        ## find left and right percentiles
        bottom,y = findIntersection(
            edges[1:],
            cumulative,
            left_percentile)

        top,y = findIntersection(
            edges[1:],
            cumulative,
            right_percentile)

        top = 10**top
        bottom = 10**bottom

        maxden = top
        dynrange = top/bottom

        if ax is not None:
            ax.step(10**edges[1:],h/h.sum()/(edges[1]-edges[0]))
            ax.text(top*1.05,0.5,'maxden',rotation=90,va='center')
            ax.plot([bottom,top],[0.5,0.5])
            ax.axvline(bottom,c='C1',ls='--',alpha=0.25)
            ax.axvline(top,c='C1')
            ax.text(np.sqrt(bottom*top),0.5/1.1,'dynrange',ha='center')
            nameAxes(ax,None,"'den' (L$_\odot$ kpc$^{-2}$)","1/N dN/d('den')",logflag=(1,0),
                supertitle="maxden=%.2g\ndynrange=%2d"%(maxden,dynrange))
            ax.get_figure().set_dpi(120)

            this_maxden = top
            if self.maxden is not None:
                #ax.text(self.maxden*1.05,0.5,'maxden',rotation=90,va='center')
                ax.axvline(self.maxden,c='C2',ls='--',alpha=0.25)
                ax.axvline(self.maxden,c='C2')
                this_maxden = self.maxden

            if self.dynrange is not None:
                ax.plot([this_maxden/self.dynrange,this_maxden],[0.4,0.4],c='C2')
                #ax.text(np.sqrt(this_maxden**2/self.dynrange),0.4/1.1,'dynrange',ha='center')



        return maxden,dynrange

    def plotParameterGrid(
        self,
        dynrange_init=100,
        maxden_init=1e-2,
        dynrange_step=None,
        maxden_step=None,
        nsteps=4,
        loud=False,
        **kwargs):
        """ Plots a grid of images that steps in maxden and dynrange to help the user decide
            what values to use, based on their aesthetic preference. Step is applied
            multiplicatively. 

            Input:

                dynrange_init = 100 -- initial value in top left corner of grid
                maxden_init = 1e-2 -- initial value in top left corner of grid
                dynrange_step = None -- step between grid thumbnails, defaults to sqrt(10)
                maxden_step = None -- step between grid thumbnails, defaults to sqrt(10)
                nsteps = 4 -- number of steps in (square) thumbnail grid
                loud = False -- flag for set_ImageParams
                **kwargs -- kwargs passed to set_ImageParams
            
            Ouput:
                
                fig - the matplotlib figure drawn to
                axs - the matplotlib axes in the grid

Example usage:
```python
starStudio.plotParameterGrid(
    dynrange_init=1e3,
    maxden_init=1,
    dynrange_step=.1,
    maxden_step=.1,
    nsteps=2,
    use_colorscheme_nasa=False)
```"""
        
        ## initialize the steps for each thumbnail in the grid
        dynrange_step = 1/np.sqrt(10) if dynrange_step is None else dynrange_step
        maxden_step = 1/np.sqrt(10) if maxden_step is None else maxden_step

        ## initialize the figure and axes
        fig,axs = plt.subplots(nrows=nsteps,ncols=nsteps)

        ## loop through the grid and set the parameters according to the initial
        ##  parameters and their steps
        for i in range(axs.shape[0]):
            for j in range(axs.shape[1]):
                ax = axs[i,j]

                ## compute this step's parameters
                dynrange = dynrange_init * dynrange_step**i
                maxden = maxden_init * maxden_step**j

                ## set parameters
                self.set_ImageParams(
                    dynrange=dynrange,
                    maxden=maxden,
                    loud=True,
                    **kwargs)

                ## actual call to render
                pixels = self.render(ax)

                ## annotate the parameters on top of the thumbnail
                nameAxes(
                    ax,None,None,None,
                    supertitle='maxden=$%1gx%.2g^{%d}$\ndynrange=$%1gx%.2g^{%d}$'%(
                        maxden_init,maxden_step,j,
                        dynrange_init,dynrange_step,i),
                    subtextkwargs={'color':'white'},subfontsize=20)

        fig.set_size_inches(4*nsteps,4*nsteps)
        fig.subplots_adjust(wspace=0,hspace=0,left=0,right=1,bottom=0,top=1)
        return fig,axs
        #fig.savefig('../src/hubble_grid.pdf',pad_inches=0,bbox_inches='tight')

append_function_docstring(StarStudio,StarStudio.set_ImageParams)
append_function_docstring(StarStudio,StarStudio.get_mockHubbleImage)
append_function_docstring(StarStudio,StarStudio.render)
append_function_docstring(StarStudio,StarStudio.predictParameters)
append_function_docstring(StarStudio,StarStudio.plotParameterGrid)
append_function_docstring(StarStudio,Studio)

##### Image projection stuff
## Stellar light attenuation projection
def raytrace_ugr_attenuation(
    x,y,z,
    mstar,ages,metals,
    h_star, 
    gx,gy,gz,
    mgas, gas_metals,
    h_gas,
    kappas,lums,
    pixels = 1200,
    xlim = None, ylim = None, zlim = None,
    QUIET=False,
    ):

    ## setup boundaries to cut-out gas particles that lay outside
    ## range
    if xlim is None:
        xlim = [np.min(x),np.max(x)]
    if ylim is None:
        ylim = [np.min(y),np.max(y)]
    if zlim is None:
        zlim = [np.min(z),np.max(z)]

    return stellar_raytrace(
        x,y,z,
        h_star,
        gx,gy,gz,
        mgas,gas_metals,
        h_gas,
        kappas,lums,
        xlim=xlim,ylim=ylim,zlim=zlim,
        pixels=pixels,
        QUIET=QUIET) 

__doc__  = ''
__doc__ = append_string_docstring(__doc__,StarStudio)
