import numpy as np 
import os,sys,getopt
import itertools

from firestudio.studios.star_studio import StarStudio


def renderStarGalaxy(
    ax,
    snapdir,snapnum,
    datadir,
    frame_half_width,frame_depth,
    edgeon=False,
    **kwargs):
    """
    Input:
        snapdir - location that the snapshots live in
        snapnum - snapshot number

        frame_half_width - half-width of image in data space
        frame_depth - half-depth of image in data space 

------- StarStudio
    Optional:
        overwrite=0 - flag to overwrite the intermediate grid files
        datadir=None - directory to output the the intermediate grid files and output png

        maxden=1e-3 - controls the saturation of the image in a non-obvious way
        dynrange=1.001 - controls the saturation of the image in a non-obvious way
        color_scheme_nasa=True - flag to use nasa colors (vs. SDSS if false) 

        single_image=None - string, if it's "Density" it will plot a column 
            density projection, if it's anything else it will be a mass weighted
            `quantity_name` projection. None will be a "two-colour" projection
            with hue determined by `quantity_name` and saturation by density

        quantity_name='Temperature' - the name of the quantity that you're mass weighting
            should match whatever array you're passing in as quantity

        snapdict=None - Dictionary-like holding gas snapshot data, open from disk if None
        star_snapdict=None - Dictionary-like holding star snapshot data, open from disk if None
        intermediate_file_name="proj_maps" ##  the name of the file to save maps to

------- Studio
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
    """
    ## make the call to renderWrapper
    render(snapdir,snapnum,edgeon,kwargs,ax)


######## Meat and Potatoes render looping functions
def renderWrapper(args):
    return render(*args)

def render(
    snapdir,snapnum,
    datadir,
    frame_half_width,frame_depth,
    edgeon,
    kwargs,
    ax=None):

    image_names = ['out_u','oul_g','out_r']
    starStudio = StarStudio(
        snapdir=snapdir,
        snapnum=snapnum,
        datadir=datadir,
        frame_half_width=frame_half_width,
        frame_depth=frame_depth,
        **kwargs)

    if edgeon:
        starStudio.renderFaceAppendEdgeViews(image_names)
    else:
        if ax is None:
            ax = plt.gca()
            ## set figure size to square
            ax.get_figure().set_size_inches(6,6)
        starStudio.render(ax,image_names)
        if ax is None:
            plt.clf()
    
def main(
    snapdir,
    snapstart,snapmax,
    datadir,
    frame_half_width,
    frame_depth,
    edgeon=0,
    **kwargs):

    if 'multiproc' in kwargs and kwargs['multiproc']:
        ## map a wrapper to a pool of processes
        argss = itertools.izip(
            itertools.repeat(snapdir),
            range(snapstart,snapmax+1),
            itertools.repeat(datadir),
            itertools.repeat(frame_half_width),
            itertools.repeat(frame_depth),
            itertools.repeat(edgeon),
            itertools.repeat(kwargs)
            itertools.repeat(None))
        my_pool = multiprocessing.Pool(int(kwargs['multiproc']))
        my_pool.map(renderWrapper,argss)
    else:
        ## just do a for loop
        for snapnum in range(snapstart,snapmax+1):
            render(snapdir,snapnum,edgeon,kwargs,ax)

if __name__=='__main__':
    argv = sys.argv[1:]
    opts,args = getopt.getopt(argv,'',[
        'snapdir=',
        'snapstart=','snapmax=',
        'pixels=','frame_half_width=','frame_depth=',
        'theta=','phi=','psi=',
        'edgeon=',
        'dynrange=','maxden=',
        'datadir=',
        'noaxis=',
        'multiproc=',
        'extract_galaxy=',
        'ahf_path=',
        'figure_label=',
        'overwrite=',
    ])    #options:
    # -r/s = use readsnap or use single snapshot loader
    #--snapdir: place where snapshots live
    #--snapstart : which snapshot to start the loop at
    #--snapmax : which snapshot to end the loop at
    #--frame_half_width : half width of frame in kpc
    #--frame_depth : half depth of frame in kpc
    #--datadir: place to output frames to

    #--theta,phi,psi : euler angles for rotation

    #--dynrange : TODO unknown
    #--maxden : TODO unknown

    #--edgeon : flag for sticking a 90 degree edge on rotation underneath 
    #--pixels : how many pixels in each direction to use, defaults to 1200
    #--noaxis : flag for removing axis and whitespace for just the pretty part
    #--multiproc : how many processes should be run simultaneously, keep in mind memory constraints
    #--extract_galaxy=False : flag to use abg_python.cosmoExtractor to extract main halo
    #--ahf_path : path relative to snapdir where the halo files are stored
    #--figure_label: text to put in the upper right hand corner

    #--overwrite: flag to  overwrite the cached projection if it exists


    for i,opt in enumerate(opts):
        if opt[1]=='':
            opts[i]=('mode',opt[0].replace('-'))
        else:
            key = opt[0].replace('-','')
            if key in ['snapdir','datadir']:
                value= opt[1]
            else:
                # turn arguments from strings to whatever
                try:
                    value = eval(opt[1])
                except:
                    value = opt[1]
            opts[i]=(key,value)
    main(**dict(opts))
