from setuptools import Extension
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup( name='PyKinematicalBroadening',
       version='0.0.8',
       description='Extragalactic Kinematics is an exciting tool that utilizes a kernel (e.g., Gaussian) to broaden models in velocity space, resulting in a highly accurate and detailed output. With this repository, you can easily apply kinematical broadening to your models and gain valuable insights into extragalactic kinematics.',
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
    
