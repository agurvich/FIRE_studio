from abg_python.plot_utils import plt

from .studio import Studio

class Composition(object):

    saveFigure = Studio.saveFigure

    def __repr__(self):
        return "Composition of (%s)"%",".join(
            [which_studio.__name__ for 
            which_studio in self.studios])

    def __init__(
        self,
        studios_tuple:tuple,
        *args,
        savefig:str=None,
        size_inches:tuple=None,
        ncols:int=2,
        nrows:int=None,
        subplots_kwargs:dict=None,
        studio_kwargss:dict=None,
        **kwargs):
        """ Make a single matplotlib figure with subplots contianing different studios with the same view.

        :param studios_tuple: tuple of :class:`firestudio.studios.studio.Studio` sub-classes
        :type studios_tuple: tuple
        :param savefig: name of file that final image should be written to. If ``None`` no file is written, defaults to ``None``
        :type savefig: str, optional
        :param size_inches: tuple of ``(width,height)`` size of figure in inches, defaults to ``None``
        :type size_inches: tuple, optional
        :param ncols: number of columns to plot side-by-side, defaults to ``2``
        :type ncols: int, optional
        :param nrows: 
            number of rows, if ``None`` calculates as many rows as necessary given\
            ``len(studios_tuple)`` and ``ncols`` defaults to ``None``
        :type nrows: int, optional
        :param subplots_kwargs: keyword arguments to be passed to ``fig.subplots_adjust``, defaults to ``None``
        :type subplots_kwargs: dict, optional
        :param studio_kwargss: 
            list of keyword argument dictionaries to be passed to each studio\
            in ``studios_tuple`` (in order) at initialization, defaults to ``None``
        :type studio_kwargss: dict, optional
        :raises ValueError: if ``len(studio_kwargss) != len(studios_tuple)``
        :raises ValueError: if ``nrows*ncols < len(studios_tuple)``
        :raises TypeError: if ``studios_tuple`` is not a tuple of length ``>=2``
        """

        nstudios = len(studios_tuple)
        if studio_kwargss is None: studio_kwargss = [{} for i in range(nstudios)]
        elif len(studio_kwargss) != nstudios: raise ValueError(
            "studio_kwargss must be the same length as studios_tuple "+
            "(%d,%d)"%(len(studio_kwargss),nstudios))
        
        for i,studio_kwargs in enumerate(studio_kwargss): 
            ## overwrite whatever savefig settings they had, 
            ##  we don't want to output individual frames
            studio_kwargs['savefig'] = None
            ## only put the scalebar in the left-most panel
            if i!= 0: studio_kwargs['scale_bar'] = False
            ## only put the label in the right-most panel
            if i!=(nstudios-1): studio_kwargs['figure_label']=''

        if ncols is None: ncols = 2
        if nrows is None: nrows= nstudios//ncols+nstudios%ncols

        if  nrows*ncols < nstudios:
            raise ValueError(
                "nrows*ncols must be >= nstudios "+
                "[(%d,%d) %d]"%(nrows,ncols,nstudios))

        self.nrows = nrows
        self.ncols = ncols

        self.size_inches = size_inches if size_inches is not None else (4*ncols,4*nrows)

        if len(studios_tuple) < 2: 
            raise TypeError(
                "studios_tuple must be a tuple of Studio"
                " sub-classes, not %s"%repr(studios_tuple))

        ## overwrite any passed kwargs with individual kwargs
        ##  for each studio specified by the studio_kwargss argument
        self.studios = [
            which_studio(*args,**{**kwargs,**which_kwargs}) for 
            which_studio,which_kwargs in 
            zip(studios_tuple,studio_kwargss)]

        self.subplots_kwargs = subplots_kwargs
        self.savefig = savefig

        self.noaxis = True
        for which_studio in self.studios: 
            self.noaxis = self.noaxis and which_studio.noaxis

        self.datadir = self.studios[0].datadir
    
    def set_ImageParams(self,**kwargs):
        """ Calls each studio's ``set_ImageParams`` method to keep them synchronized. """
        for which_studio in self.studios:
            which_studio.set_ImageParams(**kwargs)
    
    def render(self,axs=None,**kwargs):
        """ Renders each studio to each axis.

        :param axs: 
            list of matplotlib axes to render to, if ``None``\
            creates a new matplotlib figure, defaults to ``None``
        :type axs: list, optional
        :return: 
            |  ``axs`` - list of matplotlib subplot axes
            |  ``ims`` - list of image arrays
        :rtype: list, list
        """

        if axs is None: fig,axs = plt.subplots( nrows=self.nrows,ncols=self.ncols)

        fig.subplots_adjust(**self.subplots_kwargs)
        fig.set_size_inches(self.size_inches)
        
        ims = [which_studio.render(ax,**kwargs)[1]
            for ax,which_studio in zip(axs.flatten(),self.studios)]

        ## save the image
        if self.savefig is not None: self.saveFigure(fig,self.savefig)

        return axs,ims