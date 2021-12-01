from ..utils.camera_utils import Camera
try: 
    from firefly.data_reader import Settings
    firefly_enabled = True
except: 
    print("Can't import firefly, install at alexbgurvi.ch/Firefly for full features.")
    firefly_enabled = False

def camera_from_settings(filepath):
    if firefly_enabled:
        settings = Settings()
        settings.loadFromJSON(filepath,loud=False)
        ## TODO quaternion seems to be rotated by 90 degrees, weird.
        quaternion = settings['quaternion']
        print('quaternion might be 90 degrees off because webgl thinks x=z?')
        camera_focus = settings['center']
        camera_pos = settings['camera']
    else: raise NotImplementedError("Download firefly or implement a simple .json reader here and submit a PR! :)")

    return Camera(camera_pos,camera_focus,None,quaternion)