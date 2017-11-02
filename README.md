# FIRE_studio
Movie Making Utilities for FIRE simulations

This git repository was lovingly made.

## Importing to a preexisting script
After adding FIRE_studio to your python path, with a simple 
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
