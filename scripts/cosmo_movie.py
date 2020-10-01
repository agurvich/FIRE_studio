import time
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
import gc

from abg_python.galaxy.gal_utils import Galaxy
from abg_python.all_utils import filterDictionary

from firestudio.studios.gas_studio import GasStudio
snapdir = "/home/agurvich/scratch/snaps/metal_diffusion/m12i_res7100/output"

def main():
    
    ## get orientation from main galaxy
    snapnum = 600 
    galaxy = Galaxy(
        'm12i_res7100',
        snapdir,
        snapnum,
        datadir='/home/agurvich/scratch/data/metal_diffusion')

    try:
        phi_TB = galaxy.metadata.star_extract_phi_TB
        theta_TB = galaxy.metadata.star_extract_theta_TB
    except AttributeError:
        galaxy.extractMainHalo(save_meta=True)


    print(phi_TB,theta_TB)

    ## determine which snapshots we should be running on
    ##  (all available ones)
    snaps = []
    for fname in os.listdir(galaxy.snapdir):
        snaps += [int(fname.split('snapdir_')[1])]

    snaps.sort()
    snaps = np.array(snaps)

    ## check that we have all consecutive snapshots
    ##  otherwise our movie will look bad
    consecutive = np.arange(snaps[0],snaps[-1]+1)
    if (not consecutive.size == snaps.size or
        not np.all(snaps == consecutive)):
        raise IOError(snaps,consecutive)

    print(snaps)

    ## update references and clear snapshot memory
    ##  before spawning new processes
    del galaxy
    locals()
    globals()
    gc.collect()

    ## TODO wrap in multiprocessing

    with open('timing.txt','w') as handle:
        init_time = time.time()
        for snapnum in snaps[-1:]:
            makeMovieFrame(snapnum,phi_TB,theta_TB)
            handle.write("%d \t %.5f"%(snapnum,time.time()-init_time))

def makeMovieFrame(
    snapnum,
    phi_TB,
    theta_TB):

    galaxy = Galaxy(
        'm12i_res7100',
        snapdir,
        snapnum,
        datadir='/home/agurvich/scratch/data/metal_diffusion')

    galaxy.extractMainHalo(
        save_meta=True,
        force_phi_TB=phi_TB,
        force_theta_TB=theta_TB)

    fig = plt.figure()
    fig.subplots_adjust(wspace=0,hspace=0,left=0,right=1,bottom=0,top=1)

    galaxy.render(plt.gca(),loud=False)

    fig.set_size_inches(4,4)
    fig.savefig(
        'frame_%d.png'%snapnum,
        pad_inches=0,
        bbox_inches='tight',
        dpi=120)

if __name__ == '__main__':
    main()
