# FIRE_studio
Movie Making Utilities for FIRE simulations

This git repository was lovingly made.

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
`from movie_maker import renderGalaxy`, you can add a rendered galaxy to a matplotlib axis object. 
The appropriate function call resembles:
```python
renderGalaxy(
    ax,
    snapdir,snapnum,
    frame_width=gal_radius, ## half-width of frame in code units
    frame_depth=gal_radius, ## half-depth of the frame in code units
    frame_center=np.zeros(3), ## position to put the center of the frame on
    extract_galaxy=False, ## flag for whether FIRE_studio should find the main galaxy and extract it
    snapdict=snapdict, ## dictionary with snapshot info using default snapshot keys
    datadir=datadir, ## directory to save intermediate files to
    savefig=savefig, ## flag to save the image to datadir
    noaxis=noaxis, ## flag for turning off axis (1=off 0=on)
    **kwargs)
```

Where `snapdict` is a python dictionary holding the snapshot arrays with keys that match the FIRE defaults, 
`abg_python.snap_utils.openSnapshot` will do this for you. 

### additional keywords you can pass
* `theta/phi/psi` - euler angles to transform your view by
* `pixels` - resolution of image (pixels x pixels)
* `min/max_den/temp` - colorbar limits for density/temperature
* `ahf_path` - path relative to snapdir where the halo files are stored, defaults to "../ahf/halo"

## Stellar Movie Maker
---- TODO  ----
identical functionality to above exists but is not documented nor stable, see `stellar_movie_maker.py`. 
