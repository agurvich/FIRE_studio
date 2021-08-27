curdir=`pwd`
snapdir="${HOME}/snaps/metal_diffusion/m12i_res7100/output"
datadir="${HOME}/scratch/data/metal_diffusion/m12i_res7100"
snapstart=71
snapmax=600

python firestudio/gas_movie_maker.py --snapdir=${snapdir} --snapstart=${snapstart} --snapmax=${snapmax} --frame_half_width=15 --frame_half_thickness=15 --datadir=${datadir} --multiproc=20 --noaxis=1 --edgeon=1 --pixels=1200
