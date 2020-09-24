# FIRE_studio
Movie Making Utilities for FIRE simulations

This git repository was lovingly made, the code it contains is largely based off two separate (but related) visualization codes built by Volker Springel and Phil Hopkins. This is not the greatest song in the world, this is just a tribute.

requirements:
[abg_python](https://www.github.com/agurvich/abg_python)

![Sample Rendering](src/final.png)

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


gasStudio = GasStudio(
    datadir,
    snapnum,
    sim_name,
    frame_half_width=15, ## kpc
    gas_snapdict=gas_snapdict,
    )
        
gasStudio.render(
    ax,
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

image_names = ['out_u','out_g','out_r','hubble']

starStudio = StarStudio(
    datadir,
    snapnum,
    sim_name,
    frame_half_width=15, ## kpc
    gas_snapdict=gas_snapdict,
    star_snapdict=star_snapdict,
    savefig=savefig,noaxis=noaxis)
    
starStudio.render(ax,image_names)
```

Where `star_snapdict` is a python dictionary holding the snapshot arrays for `PartType4` with keys that match the FIRE defaults. `abg_python.snap_utils.openSnapshot` will do this for you. If you are making an image of an isolated galaxy, you should remember to merge the dictionaries of `PartType2`, `PartType3`, and `PartType4`.

For more information on the functionality and the different keyword arguments, see the corresponding [wiki page](https://github.com/agurvich/FIRE_studio/wiki/star_studio).
