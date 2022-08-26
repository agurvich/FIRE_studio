# FIRE_studio
[![Build Status](https://travis-ci.com/agurvich/FIRE_studio.svg?branch=main)](https://app.travis-ci.com/agurvich/FIRE_studio)
[![PyPI version](https://badge.fury.io/py/firestudio.svg)](https://badge.fury.io/py/firestudio)
<a href="https://ascl.net/2202.006"><img src="https://img.shields.io/badge/ascl-2202.006-blue.svg?colorB=262255" alt="ascl:2202.006" /></a>

Movie Making Utilities for FIRE simulations

This git repository was lovingly made, the code it contains is largely based off two separate (but related) visualization codes built by Volker Springel and Phil Hopkins. This is not the greatest song in the world, this is just a tribute.

If you use FIRE studio in a talk or paper, please acknowledge it textually and cite the code entry in the [Astrophysics Source Code Library](https://ui.adsabs.harvard.edu/abs/2022ascl.soft02006G/abstract).

Thanks, and enjoy!

## Help Policy
To receive support please fill out the issue template in the [github issues tab](https://github.com/agurvich/FIRE_studio/issues). 
Note that I only support the most recent version of FIRE studio by default and you will be required to upload a full stack trace and printed output along with the python script/pdf printout of the jupyter notebook and, if you want me to debug something I'll need pickled dictionaries with the input data as well (w/ coordinates offset, etc...) and not just a path to the snapshot hdf5 files. I ask this to help me help you!

requirements:
[abg_python](https://www.github.com/agurvich/abg_python)

![Sample Rendering](https://github.com/agurvich/FIRE_studio/blob/main/docs/source/imgs/final.png)

## Installation
To install (with ssh) clone the repository and its dependency
```bash
git clone git@github.com:agurvich/abg_python.git
git clone git@github.com:agurvich/FIRE_studio.git
```
and add them into your python path. I like to install these things to a single folder located at `${HOME}/python` that I add to my python path by including 
```bash
export PYTHONPATH="${HOME}/python:${PYTHONPATH}"
```
in my `.bashrc` file.

Then, if you'd rather install the repositories in a separate folder, you can use a soft link like so:
```bash
ln -s /path/to/repository ${HOME}/python/repository_name
```
So that you don't have to make your `PYTHONPATH` environment variable very confusing.

### Linux
You're done, congratulate yourself!

### Mac-OS / Windows
You will have to recompile the C routines in `FIRE_studio/firestudio/utils/gas_utils/HsmlAndProject_cubicSpline/` and `FIRE_studio/firestudio/utils/stellar_utils/c_libraries/`. 
My only advice is to `cd` into the directories and use `make`, if you don't know how to compile C code or end up with an error then you should focus your Google-fu efforts on intalling "Homebrew" and then use `brew install gcc` if you're on Mac-OS.
If you're on Windows then your best bet is to install Windows Subsystem for Linux (WSL) and run a virtual Linux kernel on your computer (and then you can just use `apt-get install gcc` if necessary). 

## Using FIRE_studio
There are two ways to use FIRE_studio
1) From the command line (currently broken)
2) From a Python script / Jupyter notebook

Each has its benefits/uses. If you run from within an existing Python context you can avoid having to open and reorient a snapshot (assuming you've already done that) by passing a dictionary with the required arrays. 
If you run from the command line I have included a simple multiprocessing ability so that you can render a large number of snapshots simultaneously. 


### Running from the command line
(Note that this is currently not working)
A render-loop can be started by passing any of the keyword arguments listed on the [wiki](https://github.com/agurvich/FIRE_studio/wiki) as command line arguments with the addition of `snapstart` and `snapmax`, which are the initial and final snapshots that will be rendered. 
There is also the `mps` flag that determines how many multiprocessing threads should be launched if you'd like to render the snapshots in parallel (1 thread / snapshot), make sure you have enough cores/memory for however many threads you request. 

For a gas density rendering:
`python firestudio/gas_movie_maker.py --snapdir="/home/abg6257/projects/snaps/m12i_res7100/output" --snapstart=555 --snapmax=600 --frame_width=30 --frame_depth=15 --edgeon=1 --datadir="/home/abg6257/src/FIRE_studio" --multiproc=4 --extract_galaxy=1 --noaxis=1`

or for a mock hubble (or SDSS) rendering:
`python firestudio/star_movie_maker.py --snapdir="/home/abg6257/projects/snaps/m12i_res7100/output" --snapstart=555 --snapmax=600 --frame_width=30 --frame_depth=15 --edgeon=1 --datadir="/home/abg6257/src/FIRE_studio" --multiproc=4 --extract_galaxy=1 --noaxis=1`

### Running from within a Python context
Begin by importing the studio class you would like to use, `GasStudio` for making volume renderings of the gas and its properties or `StarStudio` for mock Hubble (or SDSS) images using simulated starlight that is attenuated by dense gas (dust lanes for days!).

```python
from firestudio.studios.gas_studio import GasStudio



my_gasStudio = GasStudio(
    datadir, ## where cache files are written out to, not the simulation directory
    snapnum, ## what snapshot number, used to name cache files
    sim_name, ## what simulation name, used to name cache files
    frame_half_width=15, ## kpc, half width of image in x- and y-directions
    frame_half_thickness=15, ## kpc, half thickness of image in z-direction
    gas_snapdict=gas_snapdict, ## dictionary containing gas particle data
    star_snapdict=star_snapdict) ## dictionary containing star particle data
        
my_gasStudio.render(
    plt.gca(),
    weight_name='Masses',
    quantity_name='Temperature',
    min_weight=-0.1,
    max_weight=1.5,
    weight_adjustment_function=lambda x: np.log10(x/gasStudio.Acell)+10-6) ## msun/pc^2
    min_quantity=2,
    max_quantity=7,
    quantity_adjustment_function=np.log10)
```

Where `gas_snapdict` is a python dictionary holding the snapshot arrays for `PartType0` with keys that match the FIRE defaults. `abg_python.snap_utils.openSnapshot` will do this for you. 

For more information on the functionality and the different keyword arguments, see the corresponding [wiki page](https://github.com/agurvich/FIRE_studio/wiki/gas_studio).

```python
from firestudio.studios.star_studio import StarStudio

my_starStudio = StarStudio(
    datadir, ## where cache files are written out to, not the simulation directory
    snapnum, ## what snapshot number, used to name cache files
    sim_name, ## what simulation name, used to name cache files
    frame_half_width=15, ## kpc, half width of image in x- and y-directions
    frame_half_thickness=15, ## kpc, half thickness of image in z-direction
    gas_snapdict=gas_snapdict) ## dictionary containing gas particle data
    
my_starStudio.render(plt.gca())
```

Where `star_snapdict` is a python dictionary holding the snapshot arrays for `PartType4` with keys that match the FIRE defaults. `abg_python.snap_utils.openSnapshot` will do this for you. If you are making an image of an isolated galaxy, you should remember to merge the dictionaries of `PartType2`, `PartType3`, and `PartType4`.

For more information on the functionality and the different keyword arguments, see the corresponding [wiki page](https://github.com/agurvich/FIRE_studio/wiki/star_studio).
