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

    required_snapdict_keys = [
        'Masses',
        'Coordinates',
        'SmoothingLength',
        'Temperature',
        'AgeGyr',
        'Metallicity']

    def __repr__(self):
        return 'StarStudio instance'

####### makeOutputDirectories implementation #######
    def set_ImageParams(
        self,
        use_defaults:bool=False,
        loud:bool=True,
        **kwargs):
        """Changes the parameters of the image. 
         

        Parameters
        ----------
        use_defaults : bool, optional
            If `True` then default values of the parameters will be set 
            (potentially overwriting any previously specified parameters). 
            If `False` adjust only the keywords passed, by default False
        loud : bool, optional
            flag to print which parameters are being set/updated, by default True
            
        Keywords
        --------
        maxden : float
             controls the saturation of the image,
            sets the upper limit of the "colorbar" if None uses
            the 99 %'ile of the image surface brightness, by default None
        dynrange : float
            controls the saturation of the image,
            sets the lower limit of the "colorbar" with respect to maxden,
            if None uses the dynamic range between maxden and the 10th %'ile
            of the image surface brightness, by default None
        nodust : bool
            flag for whether dust attenuantion should be ignored, by default False
        age_max_gyr : float
            maximum age in Gyr to show stellar emission from. If None then emission from all star
            particles is considered, by default None
        
        color_scheme_nasa : bool
            flag for switching between Hubble vs. SDSS false color remapping, by default True

        Example usage:
        --------------
        ```python
        starStudio.set_ImageParams(
            maxden=0.1,
            dynrange=10,
            figure_label='Hubble')
        ```
        """

        default_kwargs = {
            'maxden':None, ## 
            'dynrange':None, ## controls the saturation of the image in a non-obvious way
            'color_scheme_nasa':True,
            'no_dust':False,
            'age_max_gyr':None} ## flag to use nasa colors (vs. SDSS if false)

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
                    if self.master_loud:
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
        if self.no_dust: self.this_setup_id+='_no_dust'
        if self.age_max_gyr is not None: self.this_setup_id+='_age_max%0.3f'%self.age_max_gyr

    append_function_docstring(set_ImageParams,Studio.set_ImageParams,prepend_string='passes `kwargs` to:\n')

    def print_ImageParams(self):
        """ Prints current image setup to console."""

        default_kwargs = {
            'maxden' : 1.0e-2, ## 
            'dynrange' : 100.0, ## controls the saturation of the image in a non-obvious way
            'color_scheme_nasa' : True, ## flag to use nasa colors (vs. SDSS if false)
            'no_dust':False,
            'age_max_gyr':None} ## flag to use nasa colors (vs. SDSS if false)

        ## print the current value, not the default value
        for arg in default_kwargs:
            print(arg,'=',getattr(self,arg))

        ## call the super class' print image params
        super().print_ImageParams()

    def quick_get_mockHubbleImage(
        self,
        lums:np.ndarray=None,
        nu_effs:list=None,
        BAND_IDS:list=None,
        **kwargs
        ):
        """Approximate the `~firestudio.studios.StarStudio.get_mockHubbleImage` routine using a 2d histogram.

        Parameters
        ----------
        lums : np.ndarray, optional
            array of luminosities for each star particle. If None then luminosities are calculated in
            the bands specified by `BAND_IDS`, by default None
        nu_effs : list, optional
            list of effective frequencies corresponding to the bands of luminosities
            provides in the `lums`. If None, uses the frequencies corresponding to the
            specified `BAND_IDS`., by default None
        BAND_IDS : list, optional
            List of indices that identify which bands to model stellar emission in. 
            If None uses Sloan UGR bands `(1,2,3)`, by default None
            ```
            lambda_eff=np.array([
                ## Bolometric (?)
                        1.e-5,  
                ## SDSS u       g       r      i      z
                        3551. , 4686. , 6165., 7481., 8931.,
                ##      U       B       V      R
                        3600. , 4400. , 5556., 6940., 
                ##      I      J       H       K
                        8700., 12150., 16540., 21790.])```

        Returns
        -------
        np.ndarray
            metal_mass_map - 2d map of metal mass used to attenuate luminosities
            along LOS in pixel, in Msun/kpc^2
        np.ndarray
            out_0 - 2d map of luminosities in first band in Lsun/kpc^2
        np.ndarray
            out_1 - 2d map of luminosities in second band in Lsun/kpc^2
        np.ndarray
            out_2 - 2d map of luminosities in third band in Lsun/kpc^2
        """

        if BAND_IDS is None:
            BAND_IDS=[1,2,3] ## used if corresponding column of lums is all 0s
       
        # apply filters, rotations, unpack snapshot data, etc...
        (kappas, lums,
            star_pos, h_star,
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
        use_metadata:bool=True,
        save_meta:bool=True,
        assert_cached:bool=False,
        loud:bool=True,
        **kwargs, 
        ):
        """Projects starlight and approximates attenuatation along line of sight
            into SDSS u, g, and r bands. 

        Parameters
        ----------
        use_metadata : bool, optional
            flag to search cache for result, by default True
        save_meta : bool, optional
            flag to cache the result, by default True
        assert_cached : bool, optional
            flag to raise an exception on a cache miss, by default False
        loud : bool, optional
            whether cache hits/misses should be announced to the console, by default True

        Returns
        -------
        np.ndarray
            metal_mass_map - 2d map of metal mass used to attenuate luminosities
            along LOS in pixel, in Msun/kpc^2
        np.ndarray
            out_0 - 2d map of luminosities in first band in Lsun/kpc^2
        np.ndarray
            out_1 - 2d map of luminosities in second band in Lsun/kpc^2
        np.ndarray
            out_2 - 2d map of luminosities in third band in Lsun/kpc^2
        """

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
        def compute_mockHubbleImage(
            self:StarStudio,
            lums:np.ndarray=None,
            nu_effs:list=None,
            BAND_IDS:list=None,
            fixed_star_hsml=0.028):

            if BAND_IDS is None:
                BAND_IDS=[1,2,3] ## used if corresponding column of lums is all 0s
            
            # apply filters, rotations, unpack snapshot data, etc...
            (kappas, lums,
                star_pos, h_star,
                gas_pos , mgas , gas_metals ,  h_gas) = self.prepareCoordinates(lums,nu_effs,BAND_IDS,fixed_star_hsml)

            ## do the actual raytracing
            gas_out,out_u,out_g,out_r = raytrace_ugr_attenuation(
                star_pos[:,0],star_pos[:,1],star_pos[:,2],
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


    def prepareCoordinates(
        self,
        lums:np.ndarray=None,
        nu_effs:list=None,
        BAND_IDS:list=None,
        fixed_star_hsml:float=None):
        """Reads snapshot data from `self.snapdict` and `self.star_snapdict` in the
        imaging volume defined by `self.Xmin`,`self.Xmax`,`self.Ymin`,`self.Ymax`,`self.Zmin`,`self.Zmax`.
        Unpacks necessary arrays and provides them to the calling imaging routine. 

        Parameters
        ----------
        lums : np.ndarray, optional
            array of luminosities for each star particle. If None then luminosities are calculated in
            the bands specified by `BAND_IDS`, by default None
        nu_effs : list, optional
            list of effective frequencies corresponding to the bands of luminosities
            provides in the `lums`. If None, uses the frequencies corresponding to the
            specified `BAND_IDS`., by default None
        BAND_IDS : list, optional
            List of indices that identify which bands to model stellar emission in. 
            If None uses Sloan UGR bands `(1,2,3)`, by default None
            ```
            lambda_eff=np.array([
                ## Bolometric (?)
                        1.e-5,  
                ## SDSS u       g       r      i      z
                        3551. , 4686. , 6165., 7481., 8931.,
                ##      U       B       V      R
                        3600. , 4400. , 5556., 6940., 
                ##      I      J       H       K
                        8700., 12150., 16540., 21790.])```

        Returns
        -------
        np.ndarray(3)
            kappas : array of effective cross sections in each band
        np.ndarray(Nstar,3)
            lums :  array of luminosities in each band
        np.ndarray(Nstar,3)
            star_pos : coordinate data for emitting stars in kpc 
        np.ndarray(Ngas)
            h_star : effective smoothing length/radius for each star particle in kpc
        np.ndarray(Ngas,3)
            gas_pos : coordinate data for attenuating gas in kpc
        np.ndarray(Ngas)
            mgas : mass data for attenuating gas in 10^10 Msun
        np.ndarray(Ngas)
            gas_metals : total metal mass fraction for each gas particle
        np.ndarray(Ngas)
            h_gas : smoothing length/radius for each gas particle in kpc
        """
        
        ## unpack the star information
        ## dont' filter star positions just yet
        star_pos = self.star_snapdict['Coordinates']

        ## cull the particles outside the frame and cast to float32
        star_pos,star_mask = self.camera.project_and_clip(star_pos)

        ages = self.star_snapdict['AgeGyr']
        if self.age_max_gyr is not None:
            star_mask = np.logical_and(star_mask,ages<self.age_max_gyr)
        ages = ages[star_mask].astype(np.float32)
        if self.master_loud: print(np.sum(star_mask),'many star particles in volume')
        
        ## try opening the stellar smoothing lengths, if we fail
        ##  let's calculate them and save them to the projection 
        ##  file

        if 'SmoothingLength' in self.star_snapdict:
            Hsml = self.star_snapdict['SmoothingLength'] 
        else:
            if fixed_star_hsml is not None:
                Hsml = np.repeat(fixed_star_hsml,star_pos.shape[0]) ## kpc
            else:
                Hsml = self.get_HSML('star')
                if Hsml.size != star_pos.shape[0]:
                    Hsml = self.get_HSML('star',use_metadata=False,save_meta=True)

        ## attempt to pass these indices along
        h_star = Hsml[star_mask].astype(np.float32)

        mstar = self.star_snapdict['Masses'][star_mask].astype(np.float32)
        metals = self.star_snapdict['Metallicity']
        if len(np.shape(metals)) > 1: metals = metals[:,0]
        metals = metals[star_mask].astype(np.float32)

        ## apply frame mask to band luminosities
        if lums is not None:
            lums = lums[:,star_mask]

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

        ## cull the particles outside the frame and cast to float32
        gas_pos,gas_mask = self.camera.project_and_clip(self.gas_snapdict['Coordinates'])

        if self.master_loud: print(np.sum(gas_mask),'many gas particles in volume')

        mgas = self.gas_snapdict['Masses'][gas_mask].astype(np.float32)
        gas_metals = self.gas_snapdict['Metallicity']
        if len(np.shape(gas_metals)) > 1: gas_metals = gas_metals[:,0]
        gas_metals = gas_metals[gas_mask].astype(np.float32)

        ## set metallicity of hot gas to 0 so there is no dust extinction
        temperatures = self.gas_snapdict['Temperature'][gas_mask]
        gas_metals[temperatures>1e5] = 0

        if "SmoothingLength" not in self.gas_snapdict:
            h_gas = self.get_HSML('gas')[gas_mask]
        else:
            h_gas = self.gas_snapdict['SmoothingLength'][gas_mask].astype(np.float32)
        
        if self.no_dust: 
            mgas = np.zeros(1)
            gas_pos = np.zeros((1,3))
            gas_metals = np.zeros(1)
            h_gas = np.array([1e-3])

        return (kappas, lums,
                star_pos, h_star,
                gas_pos , mgas , gas_metals ,  h_gas)

####### produceImage implementation #######
    def produceImage(
        self,
        quick:bool=False,
        **kwargs):
        """Generates a mock hubble image, along with any annotations/scale bars,
            using the stored image parameters.

        Parameters
        ----------
        quick : bool, optional
            flag to use a simple 2d histogram (for comparison or for quick iteration 
            as you choose image parameters), by default False

        Returns
        -------
        np.ndarray(npix_x,npix_y,3)
            RGB pixel array of the produced image
        """

        if not quick:
            gas_out,out_u,out_g,out_r = self.get_mockHubbleImage(**kwargs)
        else:
            gas_out,out_u,out_g,out_r = self.quick_get_mockHubbleImage(**kwargs)

        L_sfc_densities = np.concatenate([out_u,out_g,out_r])
        maxden_guess,dynrange_guess = self.predictParameters(L_sfc_densities=L_sfc_densities)
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
        left_percentile:float=0.1,
        right_percentile:float=0.99,
        L_sfc_densities:np.ndarray=None,
        ax:plt.Axes=None,
        quick:bool=False):
        """Guesses what the "best" values for maxden and dynrange are from
            the distribution of surface brightnesses in the current image. 
            Looks for the left_percentile and right_percentile and returns
            right_percentile and the distance between it and left_percentile
            (in log space). 
            
        Parameters
        ----------
        left_percentile : float, optional
            lower bound on image surface brightness percentile, by default 0.1
        right_percentile : float, optional
            upper bound on image surface brightness percentile, by default 0.99
        L_sfc_densities : np.ndarray, optional
            surface densities to find the percentiles of. If None uses 
            `~firestudio.studios.star_studio.StarStudio.get_mockHubbleImage` to 
            retrieve the surface densities of the projected mock hubble image, by default None
        ax : plt.Axes, optional
            axis to plot distribution of surface brightnesses
            with overlay of percentiles and such to, by default None
        quick : bool, optional
            flag to use `~firestudio.studios.star_studio.StarStudio.quick_get_mockHubbleImage`, by default False

        Returns
        -------
        float
            maxden : maximum surface brightness of the image
        float
            dynrange : distance between maximum and minimum surface brightness
                     in log space.
        """

        if (L_sfc_densities is None):
            ## read the luminosity maps
            if not quick:
                gas_out,out_u,out_g,out_r = self.get_mockHubbleImage()
            else:
                gas_out,out_u,out_g,out_r = self.quick_get_mockHubbleImage()

            L_sfc_densities = np.concatenate([out_u,out_g,out_r])

        ## concatenate the luminosity maps and take the log of the non-empty ones
        rats = np.log10(L_sfc_densities.flatten())
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
            nameAxes(
                ax,None,
                r"'den' (L$_\odot$ kpc$^{-2}$)",
                r"1/N dN/d('den')",
                logflag=(1,0),
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
        dynrange_init:float=1e2,
        maxden_init:float=1e-2,
        dynrange_step:float=None,
        maxden_step:float=None,
        nsteps:int=4,
        **kwargs):
        """ Plots a grid of images that steps in maxden and dynrange to help the user decide
            what values to use, based on their aesthetic preference. Step is applied
            multiplicatively. Accepts keyword arguments for 
            `~firestudio.studios.star_studio.StarStudio.set_ImageParams`

        Parameters
        ----------
        dynrange_init : float, optional
            initial value in top left corner of grid, by default 1e2
        maxden_init : float, optional
            initial value in top left corner of grid, by default 1e-2
        dynrange_step : float, optional
            step between grid thumbnails, by default sqrt(10), by default None
        maxden_step : float, optional
            step between grid thumbnails, by default sqrt(10), by default None
        nsteps : int, optional
            number of steps in (square) thumbnail grid, by default 4
        
        Keywords

        Returns
        -------
        plt.figure
            fig : the matplotlib figure drawn to
        list of plt.Axes
            axs : the matplotlib axes in the grid


        Example usage:
        --------------
        ```python
        starStudio.plotParameterGrid(
            dynrange_init=1e3,
            maxden_init=1,
            dynrange_step=.1,
            maxden_step=.1,
            nsteps=2,
            use_colorscheme_nasa=False)```
        """
        
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
