import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="firestudio",
    version="1.0.0",
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
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    install_requires=[            
          'abg_python'
      ],
    include_package_data=True,
)

