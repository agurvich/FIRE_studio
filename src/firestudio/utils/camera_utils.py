import numpy as np

from abg_python.math_utils import getThetasTaitBryan,construct_quaternion,q_to_rotation_matrix,q_mult,rotateQuaternion

class Camera(object):
    def __repr__(self):
        return "Camera(%s,%s) - %s "%(repr(self.camera_pos),repr(self.camera_focus),repr(self.quaternion))

    def __init__(self,camera_pos,camera_focus=None,camera_north=None,quaternion=None):
        camera_pos = np.array(camera_pos,ndmin=1)

        ## default to looking at the center
        if camera_focus is None: camera_focus = np.zeros(3)
        else: camera_focus = np.array(camera_focus)

        self.camera_focus = camera_focus
        self.camera_pos = camera_pos

        camera_normal = camera_pos-camera_focus
        self.camera_dist = np.linalg.norm(camera_normal)
        camera_normal = camera_normal/self.camera_dist

        if quaternion is None:

            ## rotate to match camera orientation-- there's a negative sign here
            ##  based on my definition of getThetasTaitBryan which was originally 
            ##  designed to align along angular momentum vector which pointed from 
            ##  the origin to the camera rather than from the camera to the origin
            theta, phi = getThetasTaitBryan(camera_normal)

            #https://yt-project.org/doc/_modules/yt/visualization/volume_rendering/off_axis_projection.html#off_axis_projection
            ##  the below is some voodoo from yt which figures out the "appropriate" camera north
            ##  somehow it nails it and gives intuitive results but i do not understand why
            if camera_north is None:
                vecs = np.identity(3)
                t = np.cross(vecs, camera_normal).sum(axis=1)
                ax = t.argmax()
                east_vector = np.cross(vecs[ax, :], camera_normal).ravel()
                camera_north = np.cross(camera_normal, east_vector).ravel()
            else:
                camera_north = np.array(camera_north)
                camera_north = camera_north / np.linalg.norm(camera_north)

            ## consider the rotation minus psi
            temp_quaternion = construct_quaternion([theta,phi],'xy')

            ## rotate the north vector into this temporary frame
            foo = rotateQuaternion(temp_quaternion,camera_north.reshape(1,3))[0]

            ## find the angle between the rotated north vector and the y-axis
            psi = 90-np.arctan2(foo[1],foo[0])*180/np.pi

            self.quaternion = construct_quaternion([theta,phi,psi],'xyz')

        else: self.quaternion = quaternion
        
        self.quat_rot_matrix = q_to_rotation_matrix(self.quaternion)
    
    def convolve_quaternion(self,new_quat):
        quaternion = q_mult(new_quat,self.quaternion)
        self.replace_quaternion(quaternion)
    
    def replace_quaternion(self,new_quat):
        self.quaternion = new_quat
        self.quat_rot_matrix = q_to_rotation_matrix(self.quaternion)

    def rotate_array(self,arr,offset=False):
        rotated_positions = np.array(np.matmul(
                self.quat_rot_matrix,
                arr.T).T,
            order='C',
            dtype=np.float32)

        if offset: rotated_center = np.array(np.matmul(
            self.quat_rot_matrix,
            self.camera_focus.reshape(3,1)).T,
            order='C',
            dtype=np.float32)
        else: rotated_center = 0

        return rotated_positions - rotated_center

    def clip(self,coords,vels):

        new_coords = self.rotate_array(coords,offset=True)
        new_vels = self.rotate_array(vels) if vels is not None else None

        ## then determine the camera distance from the camera focus
        ##  and take FOV = 45 degrees left + 45 degrees right i.e. 
        ##  xmin,xmax = -z,+z
        ##  ymin,ymax = -z,+z
        ##  where z is measured as the distance from the camera to the camera focus
        mask = np.logical_and(
            np.abs(new_coords[:,0])<self.camera_dist,
            np.abs(new_coords[:,1])<self.camera_dist)

        new_coords = new_coords[mask]
        if new_vels is not None: new_vels = new_vels[mask]

        return new_coords, new_vels

