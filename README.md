# FIRE_studio
Movie Making Utilities for FIRE simulations

This git repository was lovingly made, the code it contains is largely based off two separate (but related) visualization codes built by Volker Springel and Phil Hopkins. This is not the greatest song in the world, this is just a tribute.

requirements:
[abg_python](https://www.github.com/agurvich/abg_python)

## Features
The main function will extract the main  galaxy from a cosmological snapshot, project density and temperature maps for face- and edge-on 
orientations, combine the density/temperature to make an image whose color is the temperature and whose saturation is the density, 
and output an image to the snapshot directory provided (along with .hdf5 files containing the intermediate density/temperature maps). 

## Knobs
* Can loop over snapshots and render many without having to wrap by passing various `snapmax` and `snapstart` arguments
* Can change the aperature size and location of the image
* Can change where the images are output to (if you don't have write permissinon to snapshot directory)
* Can rotate the galaxy to various angles (with respect to face-on by default) using euler angles
* Can change the image resolution in number of pixels to a side 
* Can change the minimum/maximum limits for density/temperature colorscales

## Installation Instructions
Installation is as simple as cloning the repository and adding 
`/path/to/FIRE_studio `
to your `PYTHONPATH` environment variable. 

## Running from the command line
A render-loop can also be started with the command:

`python firestudio/gas_movie_maker.py --snapdir="/home/abg6257/projects/snaps/m12i_res7100/output" --snapstart=555 --snapmax=601 --frame_width=30 --frame_depth=15 --edgeon=1 --datadir="/home/abg6257/src/FIRE_studio" --multiproc=4 --extract_galaxy=1 --noaxis=1`


### available keywords from the command line
* `snapdir` - place where snapshots live
* `snapstart` - initial snapshot of the render-loop
* `snapmax` - final snapshot of the render-loop
* `frame_width` - half-width of frame in code units
* `frame_depth` - half-depth of frame in code units
* `datadir` - place to output images/projections to
* `edgeon` - plot 90 degree rotation underneath
* `theta/phi/psi` - euler angles to transform your view by
* `pixels` - resolution of image (pixels x pixels)
* `min/max_den/temp` - colorbar limits for density/temperature
* `multiproc` - uses that many multiprocessing threads
* `extract_galaxy` - flag to use abg_python.cosmoExtractor to extract main halo
* `ahf_path` - path relative to snapdir where the halo files are stored, defaults to "../ahf/halo"


## Running from within a python script
With a simple 
`from firestudio.studios.gas_studio import GasStudio`, you can add a rendered galaxy to a matplotlib axis object. The appropriate function call resembles:
```python
image_names = [
    'columnDensityMap',
    'massWeightedTemperatureMap',
    'two_color']

gasStudio = GasStudio(
    snapdir,snapnum,
    datadir=datadir,
    frame_half_width=radius,
    extract_galaxy=False, ## already extracted the galaxy
    snapdict=gas_snapdict,
    savefig=savefig,noaxis=noaxis,
    Hsml = gas_snapdict['SmoothingLength'])
        
gasStudio.render(ax,image_names)
```

Where `gas_snapdict` is a python dictionary holding the snapshot arrays with keys that match the FIRE defaults, `abg_python.snap_utils.openSnapshot` will do this for you. 

### additional keywords you can pass
* `theta/phi/psi` - euler angles to transform your view by
* `pixels` - resolution of image (pixels x pixels)
* `min/max_den/temp` - colorbar limits for density/temperature
* `ahf_path` - path relative to snapdir where the halo files are stored, defaults to "../ahf/halo"

## using `gas_studio.render`
### controlling two-color vs. single color
`image_names` is supposed to control this but it seems to do very little.
#### passing a single string
passing a string name will plot an image with that name and save it with that name. 
if that name is 

## using `star_studio.render`
image_names = ['out_u','out_g','out_r','hubble']


## Studio kwargs
### paths and which snapshot to open
* `snapdir`,`snapnum` - snapshot directory and snapshot number
* `datadir` - directory to put intermediate and output files
* `overwrite = False` - flag to overwrite intermediate flags
* `h5prefix=''` - string to prepend to projection file
* `this_setup_id = None` - string identifier in the intermediate projection file
* `intermediate_file_name = "proj_maps"` - the name of the file to save maps to
* `savefig = True` - save the image as a png
* `extract_galaxy = False` - uses halo center to extract region around main halo
* `ahf_path = None` - path relative to snapdir where the halo files are stored

### frame setup
* `frame_half_width`, half-width of image in x direction
* `frame_depth`, z-depth of image (thickness is 2 * frame_depth)
* `frame_center = np.zeros(3)`, center of frame in data space
* `theta=0`,`phi=0`,`psi=0`, euler rotation angles

### image parameters
* `aspect_ratio = 1` - shape of image, y/x, multiplies frame_half_width for y dimension
* `pixels = 1200` - pixels in x direction

### image annotation
* `fontsize = 12` - font size of scale bar and figure label
* `figure_label = ''` - string to be put in upper right corner
* `scale_bar = True` - flag to plot length scale bar in lower left corner
* `noaxis = True` - turns off axis ticks

## `GasStudio` kwargs
### required positional arguments passed to `Studio`
* `snapdir`,`snapnum` - snapshot directory and snapshot number
* `datadir` - directory to put intermediate and output files
* `frame_half_width` - half-width of image in x direction
* `frame_depth` - z-depth of image (thickness is 2 * frame_depth)
### color scale controls
* `min_den=-0.4` - the minimum of the density color/saturation scale (in log(n/n_units))
* `max_den=1.6` - the maximum of the density color/saturation scale (in log(n/n_units))
* `min_quantity=2` - the minimum of the quantity color scale
* `max_quantity=7` - the maximum of the quantity color scale (in log(Q/Q_units))
* `cmap='viridis'` - what colormap to use
### quantity control
* `single_image = None` - string to determine what sort of 1-color image to make
* `quantity_name='Temperature'` - quantity to make a mass weighted map/2 color image
* `take_log_of_quantity=True` - take log of mass weighted quantity map?
* `use_colorbar = False` - flag to put a colorbar
### preopened data control
* `Hsml = None` - provide smoothing lengths to speed up C routine
* `use_hsml = True` - flag to use the provided smoothing lengths (if passed)
* `snapdict = None` - provide an open snapshot dictionary to save time opening

## `StarStudio` kwargs
### required Positional Arguments
* `snapdir`,`snapnum` - snapshot directory and snapshot number
* `datadir` - directory to put intermediate and output files
* `frame_half_width` - half-width of image in x direction
* `frame_depth` - z-depth of image (thickness is 2 * frame_depth)
### color scale controls
* `maxden = 1.0e-2` - controls the saturation of the image in a non-obvious way
* `dynrange = 100.0` - controls the saturation of the image in a non-obvious way
* `color_scheme_nasa = True` - flag to use nasa colors (vs. SDSS if false)
### preopened data control
* `star_snapdict = None` - provide an open snapshot dictionary to save time opening
* `snapdict = None` - provide an open snapshot dictionary to save time opening
