import numpy as np

#### MPS LOAD BALANCING
def split_into_n_approx_equal_chunks(snap_pairs,nchunks):
    
    nrenders = len(snap_pairs)
    ## split matching pairs into groups
    indices = split_pairs(snap_pairs[:])
    splitted = np.array_split(snap_pairs,indices)
    
    ## determine how many of each matching pair there are
    n_renders_per_split = [len(this_split) for this_split in splitted]
    
    mps_chunks = []
    this_chunk = 0
    # 1 would bias early, do this if the remainder would dump a lot of pairs on the last
    ##  process, i.e. when the remainder > nchunks/2
    per_chunk = nrenders//nchunks + (nrenders%nchunks > nchunks//2)
    
    for i in range(len(indices)):
        #print(this_chunk,per_chunk)
        if this_chunk >= per_chunk:
            mps_chunks+=[indices[i-1]]
            this_chunk=0
        this_chunk+=n_renders_per_split[i]
    
    mps_chunks = np.array(mps_chunks)#-indices[0]
    mps_chunks=list(mps_chunks)
    print('split into:',np.diff([0]+mps_chunks+[len(snap_pairs)]))   
    return mps_chunks

def split_pairs(snap_pairs):
    indices = []
    changed = split_head(snap_pairs)
    cur_index = 0
    while changed > 0:
        changed = split_head(snap_pairs[cur_index:])
        cur_index+=changed
        indices+=[cur_index]

    return indices[:-1]#np.array_split(snap_pairs,indices[:-1])
    
def split_head(snap_pairs):
    split_index = np.argmax(np.logical_not(np.all(snap_pairs==snap_pairs[0],axis=1)))
    return split_index