import numpy as np

from abg_python.array_utils import filterDictionary
from abg_python.plot_utils import plt

from .studio import Studio

from matplotlib.image import imread
from tempfile import NamedTemporaryFile

class SimpleStudio(Studio):

    required_snapdict_keys = ['Coordinates']

    def produceImage(self,**kwargs):
        (xs,ys,zs) = self.__prepareCoordinates(**kwargs)

        fig,ax = plt.subplots(nrows=1,ncols=1)
        fig.subplots_adjust(wspace=0,hspace=0,left=0,right=1,bottom=0,top=1)
        fig.set_size_inches(self.npix_x/300,self.npix_y/300)

        ax.scatter(
            xs,
            ys,
            alpha=0.2,
            rasterized=True,
            s=1)

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
