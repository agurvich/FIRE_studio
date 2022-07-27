import numpy as np

from abg_python.array_utils import filterDictionary
from abg_python.plot_utils import plt

from .studio import Studio

from matplotlib.image import imread
from tempfile import NamedTemporaryFile

class SimpleStudio(Studio):

    required_snapdict_keys = ['Coordinates']
    """ these are minimum required keys for :func:`~firestudio.studios.simple_studio.SimpleStudio.render` function to run."""

    def produceImage(self,quick=True,alpha='column_density',scale_alpha=0.065,scale_radius=10.,**kwargs):
        """ Make a scatter plot rather than do any kind of projection.

        :param quick: If True do a hackup that mimics a better program, defaults to ``True``
        :type quick: bool, optional
        :param alpha: 
            Opacity of each scatter point.\
            Either 'column_density' for alpha to scale with column density,\
            or a number to apply to all alphas, defaults to ``'column_density'``
        :type alpha: str or float, optional
        :param scale_alpha: 
            What to scale ``alpha`` by. The default value is tuned for displaying ~50,000 particles, defaults to ``0.065``
        :type scale_alpha: float, optional
        :param scale_radius: 
                What to scale the radius by. The default value is tuned for displaying on ~100 kpc scales, defaults to ``10``
        :type scale_radius: float, optional

        :kwargs:
            * **snapdict_name** (`str`, `optional`) --\
                one of ``gas`` or ``star`` to identify which of ``self.gas_snapdict`` or ``self.star_snapdict`` to read data from.
            * **age_max_gyr** (`float`, `optional`) --\
                maximum age in Gyr to show stellar emission from. If ``None`` then emission from all star\
                particles is considered, defaults to ``None``

        :return: rgb image array
        :rtype: np.array
        """
        
        (xs,ys,zs) = self.prepareCoordinates(**kwargs)

        fig,ax = plt.subplots(nrows=1,ncols=1)
        fig.subplots_adjust(wspace=0,hspace=0,left=0,right=1,bottom=0,top=1)
        fig.set_size_inches(self.npix_x/300,self.npix_y/300)

        if quick:

            ax.scatter(
                xs,
                ys,
                alpha=0.2,
                rasterized=True,
                s=1)
        else:

            ## TO-DO for Alex:
            # Fix these variables up how you want.
            """
            particle_radius_data_units = ??? # Use smoothing length or sqrt of volume or whatever
            cmap = palettable.scientific.diverging.Berlin_3.mpl_colormap # Make this into what you want
            norm = matplotlib.colors.LogNorm( vmin=???, vmax=??? ) # Add in an matplotlib colors normalization for the colors
            particle_data_to_color_by = ??? # Could be temperature, for example
            particle_mass = ??? # Fill this in
            """

            width_in_data = self.Xmax - self.Xmin
            width_in_pixels = ax.get_window_extent().width
            pixels_to_points = fig.dpi / 72.
            radius = particle_radius_data_units * ( width_in_pixels / width_in_data ) * pixels_to_points * scale_radius
            s = ( radius )**2.

            # Colors
            colors = cmap( norm( particle_data_to_color_by ) )

            # Alpha
            if alpha == 'column_density':

                # Calculate alpha based on column density
                column_den = particle_mass / particle_radius_data_units**2.
                alpha_norm = matplotlib.colors.LogNorm(
                    vmin=np.nanmin( column_den ),
                    vmax=np.nanmax( column_den )
                )
                alpha = alpha_norm( column_den ) * scale_alpha

                # Remove invalid values
                alpha[alpha>1.] = 1.
                alpha[alpha<0.] = 0.
                colors[:,3] = alpha
            else:
                colors[:,3] = alpha

            # Plot itself
            ax.scatter(
                xs,
                ys,
                s = s,
                c = colors,
                edgecolors = 'none',
            )

        ax.set_xlim(self.Xmin,self.Xmax)
        ax.set_ylim(self.Ymin,self.Ymax)
        ax.axis('off')

        ## save to disk and use imread to get the rgb data, lol
        with NamedTemporaryFile(suffix='.png') as f:
            fig.savefig(
                f.name,
                bbox_inches='tight',
                dpi=300,
                facecolor='k')
            final_image = imread(f.name)
        final_image = final_image[::-1]

        plt.close(fig)
        return final_image

    def prepareCoordinates(
        self,
        snapdict_name='star',
        age_max_gyr=14,
        **extra_kwargs):
        """ filters coordinates and extracts any necessary field arrays.

        :param snapdict_name: 
            one of ``gas`` or ``star`` to identify which of ``self.gas_snapdict`` or ``self.star_snapdict`` to read data from.
        :type snapdict_name: str
        :param age_max_gyr: maximum age of star particles to filter in Gyr, defaults to ``14``
        :type age_max_gyr: int, optional
        :return:
            |  ``xs`` - x coordinates
            |  ``ys`` - y coordinates
            |  ``zs`` - z coordinates
        :rtype: np.ndarray, np.ndarray, np.ndarray
        """

        full_snapdict_name = '%s_snapdict'%snapdict_name
        ## use the masked version of the snapdict if it was passed
        if hasattr(self,'masked_'+full_snapdict_name):
            if self.master_loud:
                print("Used masked_snapdict, delete it if you don't want it anymore")
            full_snapdict_name = 'masked_'+full_snapdict_name

        snapdict = getattr(self,full_snapdict_name)

        if 'AgeGyr' in snapdict:
            snapdict = filterDictionary(snapdict,snapdict['AgeGyr']<age_max_gyr)

        return snapdict['Coordinates'].T
