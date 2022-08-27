import numpy as np

from abg_python.math_utils import rotateQuaternion_math

class Camera(object):
    def __repr__(self):
        return "Camera(%s,%s) - %s "%(repr(self.camera_pos),repr(self.camera_focus),repr(self.quaternion))

    def __init__(
        self,
        camera_pos:np.ndarray,
        camera_focus:np.ndarray=None,
        camera_up:np.ndarray=None,
        camera_roll:float=None,
        frame_half_width:float=None,
        zmin:float=0,
        zmax:float=None):
        """Class for projecting along arbitrary axes and with arbitrary focii. 

        :param camera_pos: position of camera in data-space coordinates.
        :type camera_pos: np.ndarray
        :param camera_focus: the point the camera is looking at. Along with ``camera_pos``
            defines the projection axis (``camera.camera_normal``), defaults to [0,0,0]
        :type camera_focus: np.ndarray, optional
        :param camera_up: the vector in data-space coordinates that represents "up"
            used to define the ``camera_east`` and ``camera_north`` vectors. 
            For consistent rotations, it's probably better to leave this alone and use the
            ``camera_roll`` parameter to adjust the orientation of the camera, defaults to [0,0,1]
        :type camera_up: np.ndarray, optional
        :param camera_roll: angle (in degrees) by which to rotate the camera about its normal
            vector (to spin the resulting image), defaults to None
        :type camera_roll: float, optional
        :param frame_half_width: the x-y extent of the image in the camera's frame.
            If ``None`` then the frame_half_width will be set to the distance between the camera
            position and the focus, defaults to None
        :type frame_half_width: float, optional
        :param zmin: the minimum z-distance in front of the camera a point should be for it to be rendered, 
            defaults to 0 (i.e. points behind the camera are excluded)
        :type zmin: float, optional
        :param zmax: the maximum z-distance in front of the camera before a point should be excluded,
            defaults to None (i.e. infinite)
        :type zmax: float, optional
        """

        ## for some reason this seems to work best when everything is negated?? lmfao
        camera_pos = np.array(camera_pos,ndmin=1)

        ## default to looking at the origin 
        if camera_focus is None: camera_focus = np.zeros(3)
        else: camera_focus = np.array(camera_focus)

        ## bind input
        self.camera_focus = camera_focus
        self.camera_pos = camera_pos

        ## define the look vector as pointing from the camera *to* the focus
        ##  i.e. a positive projection along camera_normal -> in front of the camera
        ## TODO why does self.camera_pos need to be negated??
        ##  camera_test() produces intuitive results only when it is :\
        self.camera_normal = self.camera_focus-(-self.camera_pos) 
        self.camera_dist = np.linalg.norm(self.camera_normal)
        self.camera_normal = self.camera_normal/self.camera_dist

        ## define frustum planes for the camera
        self.zmin = 0 if zmin is None else zmin
        ##  2x: 1x -> to get to the camera focus + another 1x to make symmetric across the camera focus
        self.zmax = 2*self.camera_dist if zmax is None else zmax
        self.frame_half_width = self.camera_dist if frame_half_width is None else frame_half_width

        ## make the the pole of our coordinate system the z-direction unless told otherwise
        if camera_up is None: camera_up = [0,0,1]

        ## find the x' and y' axes of our camera (east/north)
        self.camera_east = np.cross(camera_up, self.camera_normal).ravel()
        self.camera_east /= np.linalg.norm(self.camera_east)
        self.camera_north = np.cross(self.camera_normal, self.camera_east).ravel()

        ## roll the camera about its normal vector if requested
        ##  easiest to do with a quaternion rotation
        if camera_roll is not None:
            quat = np.zeros(4)
            quat[0] = np.cos(-camera_roll/180*np.pi/2)
            quat[1:] = np.sin(-camera_roll/180*np.pi/2)*self.camera_normal
            self.camera_east = rotateQuaternion_math(quat,self.camera_east)
            self.camera_north = rotateQuaternion_math(quat,self.camera_north)

        ## and now store the camera axes!
        self.camera_axes = np.array([
            self.camera_east,
            self.camera_north,
            self.camera_normal])        

    def project_array(self,arr:np.ndarray,offset:bool=True):
        """Apply camera projection to array, if ``offset``\
            subtract out the camera position.

        :param arr: array to project (i.e. coordinates or velocities)
        :type arr: (Npoints,3) np.ndarray
        :param offset: flag to subtract the camera position out of the array,
            this should be done for coordinate arrays but should *not* be done for
            velocities, defaults to True
        :type offset: bool, optional
        :return: projected_arr
        :rtype: np.ndarray
        """

        ## TODO why does self.camera_pos need to be negated??
        ##  camera_test() produces intuitive results only when it is :\
        if offset: offset = -self.camera_pos
        else: offset = 0

        projected_arr = np.sum(self.camera_axes*(arr-offset)[:,None],axis=-1) 

        return projected_arr


    def project_and_clip(self,coords:np.ndarray,vels:np.ndarray=None):
        """Projects both the coordinates and velocities (if provided) to the camera\
            frame and then applies the ``zmin`` and ``zmax`` clipping.

        :param coords: coordinate array to project and clip
        :type coords: np.ndarray
        :param vels: velocities of those coordinates (which should only be projected into
            the camera's frame, but not offset), defaults to None
        :type vels: np.ndarray, optional
        :return: ``new_coords``,``new_vels`` (if originally provided), ``clip_mask`` (to apply to other scalar arrays)
        :rtype: (ncoord,3) np.ndarray, [(ncoord,3) np.ndarray], (ncoord) np.ndarray
        """


        new_coords = self.project_array(coords).astype(np.float32)
        new_vels = self.project_array(vels,offset=False).astype(np.float32) if vels is not None else None

        ## then determine the camera distance from the camera focus
        ##  and take FOV = 45 degrees left + 45 degrees right i.e. 
        ##  xmin,xmax = -z,+z
        ##  ymin,ymax = -z,+z
        ##  where z is measured as the distance from the camera to the camera focus
        mask = np.logical_and(
            np.abs(new_coords[:,0])<self.frame_half_width,
            np.abs(new_coords[:,1])<self.frame_half_width)
        
        zmask = new_coords[:,-1]>self.zmin

        if self.zmax is not None: zmask = np.logical_and(zmask,new_coords[:,-1]<self.zmax)

        mask = np.logical_and(mask,zmask)

        new_coords = new_coords[mask]
        if new_vels is not None: new_vels = new_vels[mask]

        if new_vels is not None: return new_coords, new_vels, mask
        else: return new_coords, mask


def test_camera(phi=0,camera_roll=0,focus=None):

    import os
    from ..studios.studio import Studio
    from abg_python.plot_utils import plt

    offset = 0
    radius = 50
    thetas = np.linspace(1e-5,np.pi/2,12,endpoint=True)

    fig,axs = plt.subplots(nrows=3,ncols=4)
    fig.subplots_adjust(wspace=0.05,hspace=0.05)

    for ax,theta in zip(axs.flatten(),thetas[::-1]):
        camera_pos = [   radius*np.sin(theta)*np.cos(phi),
                radius*np.sin(theta)*np.sin(phi),
                radius*np.cos(theta)]


        my_camera=Camera(
            camera_pos,
            focus,
            camera_roll=camera_roll,
            frame_half_width=100)

        my_studio = Studio(
            os.path.dirname(__file__),
            0,
            'camera_test',
            camera=my_camera,
            loud=False,
            noaxis=False)

        my_studio.drawCoordinateAxes(ax,length=50)
        
        
        ax.plot(0,0,'.',c='red')
        ax.plot(-offset,-offset,'x',c='red')
        ax.set_aspect(1)

    fig.set_dpi(200)
    fig.set_size_inches(6*2,4*2)

    return fig,axs
