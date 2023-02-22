from setuptools import Extension
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup( name='PyKinematicalBroadening',
       version='0.0.5',
       description='Kinematical broadening in velocity space (km/s)',
       long_description=long_description,      # Long description read from the the readme file
       long_description_content_type="text/markdown",
       author='Jean Gomes',
       author_email='antineutrinomuon@gmail.com',
       url='https://github.com/neutrinomuon/PyKinematicalBroadening',
       install_requires=[ 'numpy','matplotlib' ],
       classifiers=[
           "Programming Language :: Python :: 3",
           "Operating System :: OS Independent",
                   ],
       package_dir={"PyKinematicalBroadening": "src/python"},
       packages=['PyKinematicalBroadening'],
       data_files=[('', ['version.txt']),],
      )
    
