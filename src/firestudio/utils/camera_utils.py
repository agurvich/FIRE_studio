import numpy as np

from abg_python.math_utils import rotateEuler,getThetasTaitBryan,construct_quaternion,q_to_rotation_matrix,q_mult

class Camera(object):
    def __init__(self,camera_pos,camera_focus=None,camera_north=None,quaternion=None):
        camera_pos = np.array(camera_pos,ndmin=1)

        ## default to looking at the center
        if camera_focus is None: camera_focus = np.zeros(3)
        else: camera_focus = np.array(camera_focus)

        if quaternion is None:
            if camera_north is not None: psi = np.arctan2(camera_north[1],camera_north[0])*180/np.pi
            else: psi = 0

            self.camera_pos = camera_pos
            self.camera_focus = camera_focus
            self.camera_north = camera_north

            ## rotate to match camera orientation-- there's a negative sign here
            ##  based on my definition of getThetasTaitBryan which was originally 
            ##  designed to align along angular momentum vector which pointed from 
            ##  the origin to the camera rather than from the camera to the origin
            theta, phi = getThetasTaitBryan(camera_pos-camera_focus)

            self.quaternion = construct_quaternion([psi,theta,phi],'zxy')
        else: self.quaternion = quaternion
        
        self.quat_rot_matrix = q_to_rotation_matrix(self.quaternion)
    
    def convolve_quaternion(self,new_quat):
        self.quaternion = q_mult(new_quat,self.quaternion)
        self.quat_rot_matrix = q_to_rotation_matrix(self.quaternion)

    def rotate_array(self,arr,offset=0):
        return np.array(np.matmul(self.quat_rot_matrix,arr.T).T,order='C') - offset

    def rotate_coordinates(self,coords):
        return self.rotate_array(coords,self.camera_focus)

    def rotate_velocities(self,vels):
        return self.rotate_array(vels)

    def clip(self,coords,vels):

        new_coords = self.rotate_coordinates(coords)
        new_vels = self.rotate_velocities(vels) if vels is not None else None

        camera_dist = np.linalg.norm(self.camera_pos-self.camera_focus)
        ## then determine the camera distance from the camera focus
        ##  and take FOV = 45 degrees left + 45 degrees right i.e. 
        ##  xmin,xmax = -z,+z
        ##  ymin,ymax = -z,+z
        ##  where z is measured as the distance from the camera to the camera focus
        mask = np.logical_and(
            np.abs(new_coords[:,0])<camera_dist,
            np.abs(new_coords[:,1])<camera_dist)

        new_coords = new_coords[mask]
        if new_vels is not None: new_vels = new_vels[mask]

        return new_coords, new_vels

