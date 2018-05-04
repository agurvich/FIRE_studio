import numpy as np 
from stellar_utils import raytrace_projection,load_stellar_hsml
import stellar_utils.make_threeband_image as makethreepic

## Stellar light attenuation projection
def calc_stellar_hsml(x,y,z):
    starhsml = load_stellar_hsml.get_particle_hsml
    h_all = starhsml(x,y,z)
    return h_all

def raytrace_ugr_attenuation(
    x,y,z,
    mstar,ages,metals,
    h_star, 
    gx,gy,gz,
    mgas, gas_metals,
    h_gas,
    pixels = 1200,
    xlim = None, ylim = None, zlim = None
    ):

    ## setup boundaries to cut-out gas particles that lay outside
    ## range
    if xlim is None:
        xlim = [np.min(x),np.max(x)]
    if ylim is None:
        ylim = [np.min(y),np.max(y)]
    if zlim is None:
        zlim = [np.min(z),np.max(z)]

    stellar_raytrace = raytrace_projection.stellar_raytrace


#   band=BAND_ID; # default=bolometric
#   j = [  0,  6,  7,  8,  9, 10, 11, 12, 13,  1,   2,   3,   4,   5] # ordering I'm used to
#   i = [  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10,  11,  12,  13] # ordering of this
#   band_standardordering = band
#   band = j[band]
#   if (band > 13): 
#       print 'BAND_ID must be < 13'; 
#       return 0;
#   
#   b=['Bolometric', \
#   'Sloan u','Sloan g','Sloan r','Sloan i','Sloan z', \
#   'Johnsons U','Johnsons B', 'Johnsons V','Johnsons R','Johnsons I', \
#   'Cousins J','Cousins H','Cousins K']
    ## pick color bands by their IDs, see above
    BAND_IDS=[9,10,11]
    #gas_out,out_u,out_g,out_r = stellar_raytrace(
    return stellar_raytrace(
        BAND_IDS,
        x,y,z,
        mstar,ages,metals,
        h_star,
        gx,gy,gz,
        mgas, gas_metals,
        h_gas,
        pixels=pixels,
        xlim=xlim,ylim=ylim,zlim=zlim,
        ADD_BASE_METALLICITY=0.001*0.02,ADD_BASE_AGE=0.0003,
        IMF_SALPETER=0,IMF_CHABRIER=1
    )
def make_threeband_image(self):


    maxden,dynrange=[(.001,1+.001),(.6,1e4)][0]
    image24, massmap = makethreepic.make_threeband_image_process_bandmaps(
        out_r,out_g,out_u,maxden=maxden,dynrange=dynrange,pixels=1200,color_scheme_nasa=1,color_scheme_sdss=0)
    plt.imshow(image24,interpolation='bicubic')#,aspect='normal')
    plt.gcf().set_size_inches(6,6)
    plt.gca().axis('off')
