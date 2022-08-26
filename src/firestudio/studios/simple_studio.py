import numpy as np

from abg_python.array_utils import filterDictionary
from abg_python.plot_utils import plt,matplotlib

from .studio import Studio

from matplotlib.image import imread
from tempfile import NamedTemporaryFile

class SimpleStudio(Studio):

    required_snapdict_keys = ['Coordinates']

    def produceImage(
        self,
        fancy=False,
        alpha='column_density',
        scale_alpha=0.065,
        scale_radius=10.,
        cmap = None,
        **kwargs):
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

        if not fancy and 'fixed_hsml' not in kwargs: kwargs['fixed_hsml'] = np.nan

        (xs,ys,zs,hs,ms,colorby) = self.__prepareCoordinates(**kwargs)
        if fancy and colorby is None: 
            if self.master_loud: print(
                "Pass in `colorby` to apply to colormap, switching fancy off.")
            fancy = False

        fig,ax = plt.subplots(nrows=1,ncols=1)
        fig.subplots_adjust(wspace=0,hspace=0,left=0,right=1,bottom=0,top=1)
        fig.set_size_inches(self.npix_x/300,self.npix_y/300)

        if not fancy:

            if colorby is None: colorby = 'C0'

            ax.scatter(
                xs,
                ys,
                alpha=0.2,
                rasterized=True,
                s=1,
                c=colorby)

        else:

            ## TO-DO for Alex:
            # Fix these variables up how you want.
            """
            norm = matplotlib.colors.LogNorm( vmin=???, vmax=??? ) # Add in an matplotlib colors normalization for the colors
            particle_data_to_color_by = ??? # Could be temperature, for example
            """

            if cmap is None: 
                import palettable
                cmap = palettable.scientific.diverging.Berlin_3.mpl_colormap # Make this into what you want

            width_in_data = self.Xmax - self.Xmin
            width_in_pixels = ax.get_window_extent().width
            pixels_to_points = fig.dpi / 72.
            radius = hs * ( width_in_pixels / width_in_data ) * pixels_to_points * scale_radius
            s = ( radius )**2.

            from matplotlib.colors import Normalize
            # Colors
            colors = cmap(Normalize(2,7,True)(np.log10(colorby)))

            # Alpha
            if alpha == 'column_density':

                # Calculate alpha based on column density
                column_den = ms / hs**2.
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
                facecolor='k',
                pad_inches=0)
            ## [::-1] because imread inherently reverses row orders
            final_image = imread(f.name)[::-1]
        final_image = np.transpose(final_image,axes=(1,0,2))

        plt.close(fig)
        return final_image

    def __prepareCoordinates(
        self,
        snapdict_name='star',
        age_max_gyr=None,
        colorby=None,
        fixed_hsml=None,
        **kwargs):

        full_snapdict_name = '%s_snapdict'%snapdict_name
        ## use the masked version of the snapdict if it was passed
        if hasattr(self,'masked_'+full_snapdict_name):
            if self.master_loud:
                print("Used masked_snapdict, delete it if you don't want it anymore")
            full_snapdict_name = 'masked_'+full_snapdict_name

        snapdict = getattr(self,full_snapdict_name)

        if age_max_gyr is not None and 'AgeGyr' in snapdict:
            snapdict = filterDictionary(snapdict,snapdict['AgeGyr']<age_max_gyr)

        coords = snapdict['Coordinates']

        ## cull the particles outside the frame and cast to float32
        coords,box_mask = self.camera.project_and_clip(coords)

        if fixed_hsml is not None: hs = np.repeat(fixed_hsml,coords.shape[0]) ## kpc
        else:
            if "SmoothingLength"  in snapdict: hs = snapdict['SmoothingLength']
            else: hs = self.get_HSML(snapdict_name)

        masses = snapdict['Masses'] if 'Masses' in snapdict.keys() else np.repeat(np.nan,coords.shape[0])

        return (
            coords[:,0],
            coords[:,1],
            coords[:,2],
            hs[box_mask],
            masses[box_mask],
            snapdict[colorby][box_mask] if colorby is not None else None)

