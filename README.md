# FIRE_studio
Movie Making Utilities for FIRE simulations

This git repository was lovingly made.

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

## Importing to a preexisting script
With a simple 
`from movie_maker import renderGalaxy`,
you can add a rendered galaxy to a matplotlib axis object. 

## Running from the command line
A render-loop can also be started with the command:
`python movie_maker.py --snapdir="/projects/b1026/agurvich/cosmo/m12i_res7000_latte/output" --snapstart=555 --snapmax=556 --frame_width=30 --frame_depth=15`


### available keywords from the command line
* `snapdir` - place where snapshots live
* `snapstart` - initial snapshot of the render-loop
* `snapmax` - final snapshot of the render-loop
* `frame_width` - half-width of frame in code units
* `frame_depth` - half-depth of frame in code units
* `datadir` - place to output images/projections to

* `theta/phi/psi` - euler angles to transform your view by
* `pixels` - resolution of image
* `min/max_den/temp` - colorbar limits for density/temperature
