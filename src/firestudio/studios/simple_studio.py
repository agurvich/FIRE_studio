import numpy as np

from abg_python.array_utils import filterDictionary
from abg_python.plot_utils import plt

from .studio import Studio

from matplotlib.image import imread
from tempfile import NamedTemporaryFile

class SimpleStudio(Studio):

    required_snapdict_keys = ['Coordinates']

    def produceImage(self,fancy=False,alpha='column_density',scale_alpha=0.065,scale_radius=10.,**kwargs):
        '''
        Args:
            fancy (bool):
                If True do a hackup that mimics a better program.

            alpha (str or float):
                Opacity of each scatter point.
                Either 'column_density' for alpha to scale with column density,
                or a number to apply to all alphas.
            
            scale_alpha (float):
                What to scale alpha by. The default value is tuned for displaying ~50,000 particles.

            scale_radius (float):
                What to scale the radius by. The default value is tuned for displaying on ~100 kpc scales.
        '''
        (xs,ys,zs) = self.__prepareCoordinates(**kwargs)

        fig,ax = plt.subplots(nrows=1,ncols=1)
        fig.subplots_adjust(wspace=0,hspace=0,left=0,right=1,bottom=0,top=1)
        fig.set_size_inches(self.npix_x/300,self.npix_y/300)

        if not fancy:

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

    def __prepareCoordinates(
        self,
        snapdict_name='star',
        age_max_gyr=14,
        **kwargs):

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
