curdir=`pwd`
python firestudio/gas_movie_maker.py --snapdir="${HOME}/snaps/m12i_res7100/output" --snapstart=600 --snapmax=600 --frame_half_width=5 --frame_depth=5 --datadir=${curdir} --multiproc=0 --extract_galaxy=1 --noaxis=1 --edgeon=0 --pixels=1200
