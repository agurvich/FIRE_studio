import numpy as np
import colorsys
import firestudio.utils.stellar_utils.utilities as util
import math
import matplotlib
matplotlib.use('Agg') ## this calls matplotlib without a X-windows GUI



# vectorize these operations for convenience below
rgb_to_hls_v = np.vectorize(colorsys.rgb_to_hls)
hls_to_rgb_v = np.vectorize(colorsys.hls_to_rgb)

def rgb_to_hls(r,g,b):
    return rgb_to_hls_v(r,g,b);

def hls_to_rgb(h,l,s):
    return hls_to_rgb_v(h,l,s);

def load_my_custom_color_tables():
    fna='heat_red'
    cdict_tmp={\
    'red':   ((0., 0.0416, 0.0416),(0.365079, 1.000000, 1.000000),(1.0, 1.0, 1.0)),\
    'green': ((0., 0., 0.),(0.71, 0.000000, 0.000000),\
            (1., 1.000000, 1.000000),(1.0, 1.0, 1.0)),\
    'blue':  ((0., 0., 0.),(0.71, 0.000000, 0.000000),(1.0, 1.0, 1.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap(fna,cdict_tmp,256);
    matplotlib.cm.register_cmap(name=fna,cmap=my_cmap);

    fna='heat_blue'
    cdict_tmp={\
    'blue':   ((0., 0.0416, 0.0416),(0.365079, 1.000000, 1.000000),(1.0, 1.0, 1.0)),\
    'green': ((0., 0., 0.),(0.365079, 0.000000, 0.000000),\
            (0.746032, 1.000000, 1.000000),(1.0, 1.0, 1.0)),\
    'red':  ((0., 0., 0.),(0.746032, 0.000000, 0.000000),(1.0, 1.0, 1.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap(fna,cdict_tmp,256);
    matplotlib.cm.register_cmap(name=fna,cmap=my_cmap);

    fna='heat_green'
    cdict_tmp={\
    'green':   ((0., 0.0416, 0.0416),(0.365079, 1.000000, 1.000000),(1.0, 1.0, 1.0)),\
    'red': ((0., 0., 0.),(0.365079, 0.000000, 0.000000),\
            (1.000000, 1.000000, 1.000000),(1.0, 1.0, 1.0)),\
    'blue':  ((0., 0., 0.),(0.1, 0.000000, 0.000000),(1.0, 1.0, 1.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap(fna,cdict_tmp,256);
    matplotlib.cm.register_cmap(name=fna,cmap=my_cmap);

    fna='heat_redyellow'
    cdict_tmp={\
    'red':   ((0., 0.0416, 0.0416),(0.365079, 1.000000, 1.000000),(1.0, 1.0, 1.0)),\
    'green': ((0., 0., 0.),(0.365079, 0.000000, 0.000000),\
            (0.746032, 1.000000, 1.000000),(1.0, 1.0, 1.0)),\
    'blue':  ((0., 0., 0.),(0.746032, 0.000000, 0.000000),(1.0, 1.0, 1.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap(fna,cdict_tmp,256);
    matplotlib.cm.register_cmap(name=fna,cmap=my_cmap);

    fna='heat_yellow'
    cdict_tmp={\
    'red':   ((0., 0.0416, 0.0416),(0.365079, 1.000000, 1.000000),(1.0, 1.0, 1.0)),\
    'green':   ((0., 0.0416, 0.0416),(0.365079, 1.000000, 1.000000),(1.0, 1.0, 1.0)),\
    'blue':  ((0., 0., 0.),(0.746032, 0.000000, 0.000000),(1.0, 1.0, 1.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap(fna,cdict_tmp,256);
    matplotlib.cm.register_cmap(name=fna,cmap=my_cmap);

    fna='heat_purple'
    cdict_tmp={\
    'red':   ((0., 0.0416, 0.0416),(0.565079, 1.000000, 1.000000),(1.0, 1.0, 1.0)),\
    'green': ((0., 0., 0.),(0.565079, 0.000000, 0.000000),\
            (0.946032, 1.000000, 1.000000),(1.0, 1.0, 1.0)),\
    'blue':   ((0., 0.0416, 0.0416),(0.565079, 1.000000, 1.000000),(1.0, 1.0, 1.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap(fna,cdict_tmp,256);
    matplotlib.cm.register_cmap(name=fna,cmap=my_cmap);

    fna='heat_orange'
    cdict_tmp={\
    'red':   ((0., 0.0416, 0.0416),(0.365079, 1.000000, 1.000000),(1.0, 1.0, 1.0)),\
    'green':   ((0., 0.0416, 0.0416),(0.365079, 0.400000, 0.400000),(1.0, 1.0, 1.0)),\
    'blue':  ((0., 0., 0.),(0.746032, 0.000000, 0.000000),(1.0, 1.0, 1.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap(fna,cdict_tmp,256);
    matplotlib.cm.register_cmap(name=fna,cmap=my_cmap);

    fna='rainbow'
    cdict_tmp = {'red':   ((0.0, 80/256., 80/256.),
                       (0.2, 0.0, 0.0),
                       (0.4, 0.0, 0.0),
                       (0.6, 256/256., 256/256.),
                       (0.95, 256/256., 256/256.),
                       (1.0, 150/256., 150/256.)),
             'green': ((0.0, 0/256., 0/256.),
                       (0.2, 0/256., 0/256.),
                       (0.4, 130/256., 130/256.),
                       (0.6, 256/256., 256/256.),
                       (1.0, 0.0, 0.0)),
             'blue':  ((0.0, 80/256., 80/256.),
                       (0.2, 220/256., 220/256.),
                       (0.4, 0.0, 0.0),
                       (0.6, 20/256., 20/256.),
                       (1.0, 0.0, 0.0))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap(fna,cdict_tmp,256);
    matplotlib.cm.register_cmap(name=fna,cmap=my_cmap);


def hue_minus_convert(hue,lightness,saturation):
    h=hue-90.; h[h<0.]+=360.; h[h>360.]-=360.;
    s=saturation*1.0; s[s<0.]=0.; s[s>1.]=1.;
    l=lightness*1.0; l[l<0.]=0.; l[l>1.]=1.;
    l=l**0.1; l=(l-np.min(l))/(np.max(l)-np.min(l));
    return h,l,s;


def invertreverse_rgb(R,G,B):
    ## invert+huereverse
    r=1.-R; g=1.-G; b=1.-B ## rgb inversion -- *different* from photoshop (cmyk inversion)
    ## convert to hls
    h,l,s=rgb_to_hls(r,g,b);
    ## huereverse
    h+=0.5; h[h>1.]-=1.; h[h<0.]+=1.;
    s*=1.1; s[s<0.]=0.; s[s>1.]=1.;
    l*=1.15; l[l<0.]=0.; l[l>1.]=1.;
    l=l**0.43; l=(l-np.min(l))/(np.max(l)-np.min(l));
    ## convert back to rgb
    r,g,b=hls_to_rgb(h,l,s);
    return r,g,b;


def temperature_map_color_index(mass_pic, temp, set_temp_max=0, set_temp_min=0, 
        huem100=0, invertreverse=0):
        
    ## custom indexing of colors for temperature maps
    ## Assuming Temperature is the weighting :: interesting temp scale is from 1.0d4 to 1.0d7
    temp_max=np.max(temp); temp_min=np.median(temp)/(np.max(temp)/np.median(temp));
    if (set_temp_max != 0): temp_max=set_temp_max;
    if (set_temp_min != 0): temp_min=set_temp_min;
    tmp_scale=np.log10(temp/temp_min)/np.log10(temp_max/temp_min);
    print('Temp. Scale Min/Max = ',np.min(tmp_scale),temp_min,np.max(tmp_scale),temp_max)
    
    N0=temp.shape[0]; image24=np.zeros([N0,N0,3]);
	## different kernels we can use for the temp weighting
	## use a simple HLS (hue, lightness, saturation) model and use that to get the RGB colors
	## hue runs 0=red - 60=yellow - 120=green - 240=blue - 320=magenta
		
    ## get something like Volker's scheme with 
    tmp_scale[tmp_scale<0.]=0.; tmp_scale[tmp_scale>1.]=1.;
    hue = 240.+tmp_scale*180.; hue[hue >= 360.]-=360.; hue[hue<0.]+=360.;
    saturation = (0.*tmp_scale + 0.99) * 0.5  ## richer colors -> 1, more gray -> 0
    L = mass_pic / 256.;
    lightness = 0.01+0.985*L;  ## 0=black, 1=white
    lightness=0.5 + (lightness-0.5)*0.75;
    lightness[lightness<0.]=0.; lightness[lightness>1.]=1.;
    if (huem100==1): hue,lightness,saturation = hue_minus_convert(hue,lightness,saturation);
    hue /= 360.; ## python uses units of 0-1 for colors
    
    R,G,B = hls_to_rgb(hue,lightness,saturation);
    if (invertreverse==1): R,G,B = invertreverse_rgb(R,G,B);

    image24[:,:,0]=R; image24[:,:,1]=G; image24[:,:,2]=B
    return image24;


class CustomLightSource(object):
    """
    Create a light source coming from the specified azimuth and elevation.
    Angles are in degrees, with the azimuth measured
    clockwise from north and elevation up from the zero plane of the surface.
    The :meth:`shade` is used to produce rgb values for a shaded relief image
    given a data array.
    (Modified version of matplotlib 'lightsource' that allows for fixed limits 
     on shading, important for consistency between images in a series)
    """
    def __init__(self,azdeg=315,altdeg=45,\
                 hsv_min_val=0,hsv_max_val=1,hsv_min_sat=1,hsv_max_sat=0):
       """
       Specify the azimuth (measured clockwise from south) and altitude
       (measured up from the plane of the surface) of the light source
       in degrees.

       The color of the resulting image will be darkened
       by moving the (s,v) values (in hsv colorspace) toward
       (hsv_min_sat, hsv_min_val) in the shaded regions, or
       lightened by sliding (s,v) toward
       (hsv_max_sat hsv_max_val) in regions that are illuminated.
       The default extremes are chose so that completely shaded points
       are nearly black (s = 1, v = 0) and completely illuminated points
       are nearly white (s = 0, v = 1).
       """
       self.azdeg = azdeg
       self.altdeg = altdeg
       self.hsv_min_val = hsv_min_val
       self.hsv_max_val = hsv_max_val
       self.hsv_min_sat = hsv_min_sat
       self.hsv_max_sat = hsv_max_sat

    def shade(self,data,cmap,vmin=None,vmax=None,limit_elevation=False):
        """
        Take the input data array, convert to HSV values in the
        given colormap, then adjust those color values
        to given the impression of a shaded relief map with a
        specified light source.
        RGBA values are returned, which can then be used to
        plot the shaded image with imshow.

        Parameters
        ----------
        data            Input data array
        vmin            Minimum data value for colormap.  Default min(data).
        vmax            Maximum data value for colormap.  Default max(data).
        limit_elevation Limit the elevation in the shading routine? Default False.
                        If true, the elevation will be limited by vmin and vmax.

        Returns 
        -------
        rgb             Shaded RGBA values, suitable for use with imshow.

        """

        if (vmin is not None) or (vmax is not None):
            limitschanged = True
        else:
            limitschanged = False
        if vmin is None:
            vmin = np.min(data)
        if vmax is None:
            vmax = np.max(data)
        rgb0 = cmap((data-vmin)/(np.float(vmax-vmin)))
        #avoid using extra memory if copy of array not needed
        if limitschanged and limit_elevation:
            d = data.copy()
            d[d<vmin] = vmin
            d[d>vmax] = vmax
            rgb1 = self.shade_rgb(rgb0, elevation=d)
        else:
            rgb1 = self.shade_rgb(rgb0, elevation=data)
        rgb0[:,:,0:3] = rgb1
        return rgb0

    def shade_rgb(self,rgb, elevation, fraction=1., vmin=None,vmax=None):
        """
        Take the input RGB array (ny*nx*3) adjust their color values
        to given the impression of a shaded relief map with a
        specified light source using the elevation (ny*nx).
        A new RGB array ((ny*nx*3)) is returned.
        """
        # imagine an artificial sun placed at infinity in
        # some azimuth and elevation position illuminating our surface. The parts of
        # the surface that slope toward the sun should brighten while those sides
        # facing away should become darker.
        # convert alt, az to radians
        az = self.azdeg*np.pi/180.0
        alt = self.altdeg*np.pi/180.0
        # gradient in x and y directions
        dx, dy = np.gradient(elevation)
        
        # clip according to the imposed limits:
        if vmin is None: vmin=np.min(dx);
        if vmax is None: vmax=np.max(dx);
        dx[dx<vmin]=vmin; dx[dx>vmax]=vmax;
        if vmin is None: vmin=np.min(dy);
        if vmax is None: vmax=np.max(dy);
        dy[dy<vmin]=vmin; dy[dy>vmax]=vmax;

        slope = 0.5*np.pi - np.arctan(np.hypot(dx, dy))
        aspect = np.arctan2(dx, dy)
        intensity = np.sin(alt)*np.sin(slope) + np.cos(alt)*np.cos(slope)*np.cos(-az -\
                aspect - 0.5*np.pi)
        # rescale to interval -1,1
        # +1 means maximum sun exposure and -1 means complete shade.
        imin = np.min(intensity)
        imax = np.max(intensity)
        imin1 = (np.sin(alt) + np.cos(alt)*vmin) / np.sqrt(1.+vmin*vmin)
        imin2 = (np.sin(alt) - np.cos(alt)*vmin) / np.sqrt(1.+vmin*vmin)
        imin3 = (np.sin(alt) + np.cos(alt)*vmax) / np.sqrt(1.+vmax*vmax)
        imin4 = (np.sin(alt) - np.cos(alt)*vmax) / np.sqrt(1.+vmax*vmax)
        imin = np.min(np.array([imin1,imin2,imin3,imin4]))
        imax = np.max(np.array([imin1,imin2,imin3,imin4]))
        print(' intensity limits being used in this map == ',
            np.min(elevation),np.max(elevation),np.min(intensity),np.max(intensity),
            np.min(dx),np.max(dx),np.min(dy),np.max(dy),np.min(slope),np.max(slope),
            np.min(aspect),np.max(aspect),imin,imax,az,alt)
        
        intensity[intensity<imin] = imin
        intensity[intensity>imax] = imax
        intensity = (intensity - imin)/(imax - imin)

        intensity = (2.*intensity - 1.)*fraction
        # convert to rgb, then rgb to hsv
        hsv = matplotlib.colors.rgb_to_hsv(rgb[:,:,0:3])
        # modify hsv values to simulate illumination.
        hsv[:,:,1] = np.where(np.logical_and(np.abs(hsv[:,:,1])>1.e-10,intensity>0),\
                (1.-intensity)*hsv[:,:,1]+intensity*self.hsv_max_sat, hsv[:,:,1])
        hsv[:,:,2] = np.where(intensity > 0, (1.-intensity)*hsv[:,:,2] +\
                intensity*self.hsv_max_val, hsv[:,:,2])
        hsv[:,:,1] = np.where(np.logical_and(np.abs(hsv[:,:,1])>1.e-10,intensity<0),\
                (1.+intensity)*hsv[:,:,1]-intensity*self.hsv_min_sat, hsv[:,:,1])
        hsv[:,:,2] = np.where(intensity < 0, (1.+intensity)*hsv[:,:,2] -\
                intensity*self.hsv_min_val, hsv[:,:,2])
        hsv[:,:,1:] = np.where(hsv[:,:,1:]<0.,0,hsv[:,:,1:])
        hsv[:,:,1:] = np.where(hsv[:,:,1:]>1.,1,hsv[:,:,1:])
        # convert modified hsv back to rgb.
        return matplotlib.colors.hsv_to_rgb(hsv)
        
