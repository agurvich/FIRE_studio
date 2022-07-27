import scipy
import numpy as np

from abg_python.array_utils import findIntersection
from abg_python.galaxy.metadata_utils import metadata_cache
from abg_python.plot_utils import plt,nameAxes

from .studio import Studio
from ..utils.stellar_utils import (
    raytrace_projection_compute, ## project along the LoS and attenuate if necessary for 3 bands simultaneously
    make_threeband_image_process_bandmaps, ## take 3 bands, overlay them, and map to RGB to make final image
    layer_band_images) 

class FIREStudio(Studio):

    required_snapdict_keys = ['Masses','Coordinates','SmoothingLength','Temperature']
    """ these are minimum required keys for :func:`~firestudio.studios.FIRE_studio.FIREStudio.render` function to run."""

    def __repr__(self):
        return 'FIREStudio instance'

    def set_ImageParams(
        self,
        use_defaults=False,
        loud=True,
        **kwargs):
        """ Changes the parameters of the image. Also calls :class:`~firestudio.studios.studio.Studio`'s\
            :func:`~firestudio.studios.studio.Studio.set_ImageParams`, passing along any unmatched kwargs.

        :param use_defaults: 
            If ``True`` then default values of the parameters will be set\
            (potentially overwriting any previously specified parameters).\
            If ``False`` adjust only the keywords passed, defaults to ``False``
        :type use_defaults: bool, optional
        :param loud: 
            flag to print which parameters are being set/updated, defaults to ``True``
        :type loud: bool, optional 

        :kwargs:
            * **tcuts** (`tuple`, `optional`) -- \
                temperature cuts for each of the three RGB bands
            * **maxden** (`float`, `optional`) -- \
                controls the saturation of the image,\
                sets the upper limit of the "colorbar" if ``None`` uses\
                the 99 %'ile of the image surface brightness, defaults to ``None``
            * **dynrange** (`float`, `optional`) -- \
                controls the saturation of the image,\
                sets the lower limit of the "colorbar" with respect to maxden,\
                if ``None`` uses the dynamic range between ``maxden`` and the 10th %'ile\
                of the image surface brightness, defaults to ``None``
            
        :Example usage:
            .. code-block:: python

                fireStudio.set_ImageParams(
                    tcuts = (300,2e4,3e5)
                    figure_label='t = 13.8 Gyr')
"""

        default_kwargs = {
            'tcuts':(300,2e4,3e5),
            'maxden':2e-6,
            'dynrange':30}

        for kwarg in list(kwargs):
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
                    and kwarg != 'this_setup_id'
                    and self.master_loud
                    and loud):
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
        super().set_ImageParams(use_defaults=use_defaults,loud=loud,**kwargs)

    set_ImageParams.default_kwargs = {
            'use_colorbar':False,
            'cbar_label':'',
            'cbar_logspace':True}
    set_ImageParams.default_kwargs.update(
        Studio.set_ImageParams.default_kwargs)

    def print_ImageParams(self):
        """ Prints current image setup to console."""

        default_kwargs = {
            'tcuts':(300,4000,10000),
            }

        ## print the current value, not the default value
        for arg in default_kwargs:
            print(arg,'=',getattr(self,arg))

        ## call the super class' print image params
        super().print_ImageParams()

    def get_gasThreebandImage(
        self,
        use_metadata=True,
        save_meta=True,
        assert_cached=False,
        loud=True,
        **kwargs, 
        ):
        """ routine to use ``raytrace_projection_compute`` to make mock gas images, with three color channels for different temperature ranges

        :param use_metadata: flag for whether a cached result should be used (if it exists), defaults to ``True``
        :type use_metadata: bool, optional
        :param save_meta: flag to save the result in the cache, defaults to ``True``
        :type save_meta: bool, optional
        :param assert_cached: flag to require a cache hit and raise an exception otherwise, defaults to ``False``
        :type assert_cached: bool, optional
        :param loud: flag for whether cache hits/misses should be announced to the console, defaults to ``True``
        :type loud: bool, optional

        :kwargs:
            * **use_log_t** (`bool`, `optional`) -- \
                defaults to ``True``
            * **isosurfaces** (`bool`, `optional`) -- \
                defaults to ``False``
            * **add_temperature_weights** (`bool`, `optional`) -- \
                defaults to ``False``

        :return: 
            | ``out_all`` -\
            | ``out_cold`` -\
            | ``out_warm`` -\
            | ``out_hot`` -
        :rtype: np.ndarray, np.ndarray, np.ndarray, np.ndarray
        """

        @metadata_cache(
            self.this_setup_id,  ## hdf5 file group name
            ['threeband_allMassMap',
            'threeband_coldMassTempMap', 
            'threeband_warmMassTempMap',
            'threeband_hotMassTempMap'],
            use_metadata=use_metadata,
            save_meta=save_meta,
            assert_cached=assert_cached,
            loud=loud,
            force_from_file=True)  ## read from cache file, not attribute of object
        def compute_gasThreebandImage(
            self,
            use_log_t=True,
            isosurfaces=False,
            add_temperature_weights=False):
            # apply filters, rotations, unpack snapshot data, etc...
            (coords,hsml,mass,
            wt1,wt2,wt3,
            kappas) = self.prepareCoordinates(
                use_log_t=use_log_t,
                isosurfaces=isosurfaces,
                add_temperature_weights=add_temperature_weights)

            ## do the actual raytracing
            out_all,out_cold,out_warm,out_hot = raytrace_projection_compute(
                coords[:,0],coords[:,1],coords[:,2],
                hsml,mass,
                wt1,wt2,wt3,
                kappas[0],kappas[1],kappas[2],
                xlim=(self.Xmin, self.Xmax),
                ylim=(self.Ymin, self.Ymax),
                zlim=(self.Zmin, self.Zmax),
                pixels=self.pixels)

            return out_all,out_cold,out_warm,out_hot
        return compute_gasThreebandImage(self,**kwargs)

    def prepareCoordinates(
        self,
        use_log_t=True,
        isosurfaces=False,
        add_temperature_weights=False):
        """ many experiments here: doing gas isosurfaces with broad kernels\
            and overlaying a custom set of color tables after the fact seems\
            best. adjust opacity (kappa_units) and kernel width as needed. \
            also, can set 'dynrange=0' to automatically produce a given \
            fraction of saturated/unsaturated pixels
            
        :param use_log_t: _description_, defaults to ``True``
        :type use_log_t: bool, optional
        :param isosurfaces: _description_, defaults to ``False``
        :type isosurfaces: bool, optional
        :param add_temperature_weights: _description_, defaults to ``False``
        :type add_temperature_weights: bool, optional
        :return: _description_
        :rtype: _type_
        """

        #if KAPPA_UNITS is None: KAPPA_UNITS = 2.0885*np.array([1.1,2.0,1.5])
        #if kernel_widths is None: kernel_widths=np.array([0.8,0.3,0.6])
        KAPPA_UNITS = 2.0885*np.array([0.5,15.0,1.5])
        kernel_widths = np.array([0.8,0.3,0.6])

        full_snapdict_name = 'gas_snapdict'
        ## use the masked version of the snapdict if it was passed
        if hasattr(self,'masked_'+full_snapdict_name):
            if self.master_loud:
                print("Used masked_snapdict, delete it if you don't want it anymore")
            full_snapdict_name = 'masked_'+full_snapdict_name

        snapdict = getattr(self,full_snapdict_name) 

        ## unpack the snapshot data from the snapdict
        coords = snapdict['Coordinates'] ## kpc

        if "SmoothingLength" not in snapdict:
            hsml = self.get_HSML('gas')
            assert type(hsml) == np.ndarray
            if 'masked_' in full_snapdict_name:
                hsml = hsml[self.mask]
        else:
            hsml = snapdict['SmoothingLength'] ## kpc

        ## rotate by euler angles if necessary
        coords = self.camera.rotate_array(coords,offset=True)

        ## cull the particles outside the frame and cast to float32
        box_mask = self.cullFrameIndices(coords)

        if self.master_loud:
            print("projecting %d particles"%np.sum(box_mask))

        coords = coords[box_mask].astype(np.float32)
        hsml = hsml[box_mask].astype(np.float32)

        gas_T = snapdict['Temperature'][box_mask].astype(np.float32)

        tcuts = self.tcuts
        if(use_log_t==1):
            tcuts=np.log10(tcuts)
            gas_T=np.log10(gas_T)

        ## ascribe weights to each gas particle based on their temperature
        ##  in each band
        ## continuous smoothing with gaussians for the temperature:
        if (isosurfaces==1):
            wt1 = np.exp(-(gas_T-tcuts[0])*(gas_T-tcuts[0])/(2.*kernel_widths[0]*kernel_widths[0]))
            wt2 = np.exp(-(gas_T-tcuts[1])*(gas_T-tcuts[1])/(2.*kernel_widths[1]*kernel_widths[1]))
            wt3 = np.exp(-(gas_T-tcuts[2])*(gas_T-tcuts[2])/(2.*kernel_widths[2]*kernel_widths[2]))
        else: ## isosurfaces==0, so do total in integral ranges set by temperature_cuts  
            wt1 = 0.5*(1.0-scipy.special.erf((gas_T-tcuts[0])/(np.sqrt(2.)*kernel_widths[0])))
            wt3 = 0.5*(1.0-scipy.special.erf((tcuts[1]-gas_T)/(np.sqrt(2.)*kernel_widths[1])))
            wt2 = 1.-wt1-wt3; wt2[wt2<0.]=0.

        # weighting by sqrt(temp) makes the contributions more similar by temperature bins
        wtfn = snapdict['Masses'][box_mask].astype(np.float32) #* np.sqrt(gas_temperature/1.0e4) 
        if (add_temperature_weights==1): wtfn *= np.sqrt(1. + gas_T/1.0e4) 

        wt1*= wtfn; wt2*=wtfn; wt3*=wtfn
        kappa = 200. * (1.+np.zeros((3)))
        kappa *= KAPPA_UNITS

        return (coords,hsml,wtfn,
            wt1,wt2,wt3,
            kappa)

    def quick_get_gasThreebandImage(self,**kwargs):
        """ Approximates :func:`~firestudio.studios.FIRE_studio.FIREStudio.get_gasThreebandImage` but using 2d histograms

        :return: _description_
        :rtype: _type_
        """

        # apply filters, rotations, unpack snapshot data, etc...
        (coords,hsml,mass,
        wt1,wt2,wt3,
        kappas) = self.prepareCoordinates()

        xedges = np.linspace(self.Xmin,self.Xmax,self.pixels)
        out_all,_,_ = np.histogram2d(
            coords[:,0],
            coords[:,1],
            bins=xedges,
            weights=mass)

        out_cold,_,_ = np.histogram2d(
            coords[:,0],
            coords[:,1],
            bins=xedges,
            weights=wt1)

        out_warm,_,_ = np.histogram2d(
            coords[:,0],
            coords[:,1],
            bins=xedges,
            weights=wt2)

        out_hot,_,_ = np.histogram2d(
            coords[:,0],
            coords[:,1],
            bins=xedges,
            weights=wt3)


        return out_all,out_cold,out_warm,out_hot

####### produceImage implementation #######
    def render(self,ax:plt.Axes=None,**kwargs):
        """ Plots a projected image using the stored image parameters.

        :param ax: axis to plot image to, if ``None`` will create a new figure, defaults to ``None``
        :type ax: plt.Axes, optional
        :raises NotImplementedError: _description_
        :return:
            |  ``ax`` -- the axis the image was plotted to
            |  ``final_image`` -- output RGB pixel array
        :rtype: plt.Axes, np.ndarray(Npix_x,Npix_y,3)

        :Example usage:
            .. code-block:: python

                fireStudio.render()
        """


        if ax is None: fig,ax = plt.figure(),plt.gca()
        else: fig = ax.get_figure()

        ## remap the C output to RGB space
        final_image = self.produceImage(**kwargs)

        ## plot that RGB image and overlay scale bars/text
        self.plotImage(ax,final_image)

        ## save the image
        if self.savefig is not None: self.saveFigure(fig,self.savefig)

        return ax,final_image

    def produceImage(
        self,
        quick=False,
        **kwargs
        ):

        if not quick: out_all,out_cold,out_warm,out_hot = self.get_gasThreebandImage(**kwargs)
        else: out_all,out_cold,out_warm,out_hot = self.quick_get_gasThreebandImage(**kwargs)

        if self.maxden is None:
            maxden_guess,dynrange_guess = self.predictParameters()
            print("setting maxden to best guess %.2g"%maxden_guess)
            self.set_ImageParams(
                maxden=maxden_guess)

        if self.dynrange is None:
            maxden_guess,dynrange_guess = self.predictParameters()
            print("setting dynrange to best guess %.2g"%dynrange_guess)
            self.set_ImageParams(
                dynrange=dynrange_guess)

        image24, massmap = make_threeband_image_process_bandmaps(
            out_cold,out_warm,out_hot,
            maxden=self.maxden,
            dynrange=self.dynrange,
            QUIET=not self.master_loud)
                
                #color_scheme_nasa=nasa_colors,
                #color_scheme_sdss=sdss_colors)

        final_image = layer_band_images(image24, massmap)

        return np.transpose(final_image,axes=(1,0,2))

    def predictParameters(
        self,
        left_percentile:float=0.1,
        right_percentile:float=0.99,
        all_bands:np.ndarray=None,
        ax:plt.Axes=None,
        quick:bool=False):
        """ Guesses what the "best" values for maxden and dynrange are from\
            the distribution of surface brightnesses in the current image.\
            Looks for the left_percentile and right_percentile and return\
            right_percentile and the distance between it and left_percentil\
            (in log space).

        :param left_percentile: 
            lower bound on image surface brightness percentile, defaults to ``0.1``
        :type left_percentile: float, optional
        :param right_percentile: 
            upper bound on image surface brightness percentile, defaults to ``0.99``
        :type right_percentile: float, optional
        :param all_bands: 
            optionally plots distribution of surface brightnesses\
            (in some units...) with overlay of percentiles and such, defaults to ``None``
        :type all_bands: _type_, optional
        :param ax: _description_, defaults to None
        :type ax: _type_, optional
        :param quick: 
                flag to use a simple 2d histogram (for comparison or for quick iteration\
                as the user defines the image parameters), defaults to ``False``
        :type quick: bool, optional
        :raises NotImplementedError: if ``quick``
        :return: 
            |  ``maxden`` -- maximum surface brightness of the image
            |  ``dynrange`` -- distance between maximum and minimum surface brightness in log space. 
        :rtype: float, float
        """

        if (all_bands is None):
            ## read the luminosity maps
            if not quick: all_bands = np.concatenate(self.get_gasThreebandImage(assert_cached=True,loud=False)[1:])
            else: raise NotImplementedError("Can't generate all_bands with quick.")

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

####### plotImage implementation #######
    ## run Studio's plotImage method
    def plotImage(self,*args,**kwargs): super().plotImage(*args,**kwargs)