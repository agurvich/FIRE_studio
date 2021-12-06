import subprocess
import os

from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.egg_info import egg_info

####
#https://stackoverflow.com/questions/19569557/pip-not-picking-up-a-custom-install-cmdclass
#BEGIN CUSTOM INSTALL COMMANDS
#These classes are used to hook into setup.py's install process. Depending on the context:
#$ pip install my-package

#Can yield `setup.py install`, `setup.py egg_info`, or `setup.py develop`

def custom_command(prepend=''):
    C_routine_subdir = 'src/firestudio/utils/C_routines'
    C_routines = os.listdir(C_routine_subdir)
    for C_routine in C_routines:
        process = subprocess.Popen(
            "make",
            shell=True,
            cwd=os.path.join(C_routine_subdir,C_routine))
        process.wait()
        
class CustomInstall(install):
    def run(self):
        custom_command()
        install.run(self)


class CustomDevelop(develop):
    def run(self):
        custom_command()
        develop.run(self)


class CustomEggInfo(egg_info):
    def run(self):
        custom_command(prepend='../../')
        egg_info.run(self)

####  END CUSTOM INSTALL COMMANDS 

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="firestudio",
    version="2.0.1",
    author = 'Alex Gurvich',
    author_email = 'agurvich@u.northwestern.edu',
    description="Rendering code for FIRE simulation data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/agurvich/FIRE_studio",
    project_urls={
        "Bug Tracker": "https://github.com/agurvich/FIRE_studio/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.6",
    install_requires=[            
          'abg_python',
          'numpy',
          'scipy',
          'matplotlib'
      ],
    include_package_data=True,
    ## to compile C code
    cmdclass={
        'install': CustomInstall,
        'develop': CustomDevelop,
        'egg_info': CustomEggInfo},
)

