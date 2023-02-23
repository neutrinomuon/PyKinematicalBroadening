### PyKinematicalBroadening
email: [antineutrinomuon@gmail.com](mailto:antineutrinomuon@gmail.com), [jean@astro.up.pt](mailto:jean@astro.up.pt)

© Copyright ®

J.G. - Jean Gomes

<hr>

![python3](https://img.shields.io/pypi/pyversions/pyfluxconserving)

<hr>

#### <b>RESUME</b>

PyKinematicalBroadening is an Extragalactic Kinematics repository for applying a kernel in velocity space to models in order to obtain the respective broadened model. This is a python code that performs kinematical broadening of a spectrum by applying a kernel in velocity space to a model, and obtaining the respective broadened model. The code defines the function broadening, which performs the convolution with a Gaussian kernel. The kernel is generated using a certain number of points, which can be set with the Ni_Gauss parameter. The code then reads in a test spectrum from a file, interpolates it onto a set of equally spaced wavelength values, and then plots the original and broadened spectra for different velocity dispersions.

In detail, the GaussianConvolution function convolves a given input spectrum fluxes_o defined at wavelengths lambda_o with a Gaussian kernel of width vd_sigma and mean velocity vc0_gals. The kernel is defined with Ni_Gauss points, which should be at least as large as vd_sigma. The output spectrum is defined at wavelengths lambda_s, and is returned as fluxes_s. The fill_val parameter defines the value to use for regions outside of the original wavelength range, and verbosity controls the level of detail of console output.

The main code reads in a test spectrum from a file and interpolates it onto a set of equally spaced wavelength values. It then loops over different velocity dispersions and calls 'broadening' for each one, broadening the spectrum and plotting the results.

<hr>

#### <b>INSTALLATION</b>

You can easily install <a href=https://pypi.org/project/PyKinematicalBroadening/>PyKinematicalBroadening</a> by using pip - <a href='https://pypi.org/'>PyPI - The Python Package Index</a>:
<pre>
<code>
pip install PyKinematicalBroadening
</code>
</pre>
or by using a generated conda repository <a href='https://anaconda.org/neutrinomuon/PyKinematicalBroadening'>https://anaconda.org/neutrinomuon/PyKinematicalBroadening</a>:

[![badgetanaconda](https://anaconda.org/neutrinomuon/PyKinematicalBroadening/badges/version.svg)](https://anaconda.org/neutrinomuon/PyKinematicalBroadening/badges/version.svg)
[![badgetreleasedate](https://anaconda.org/neutrinomuon/PyKinematicalBroadening/badges/latest_release_date.svg)](https://anaconda.org/neutrinomuon/PyKinematicalBroadening/badges/latest_release_date.svg)
[![badgetplatforms](https://anaconda.org/neutrinomuon/PyKinematicalBroadening/badges/platforms.svg
)](https://anaconda.org/neutrinomuon/PyKinematicalBroadening/badges/platforms.svg)
<pre>
<code>
conda install -c neutrinomuon PyKinematicalBroadening
</code>
</pre>
OBS.: Linux, OS-X ad Windows pre-compilations available in conda.

You can also clone the repository and install by yourself in your machine:
<pre>
<code>
git clone https://github.com/neutrinomuon/PyKinematicalBroadening
python setup.py install
</code>
</pre>

<hr>

#### <b>EXAMPLE</b>

Example of the test_spectrum test_spectrum.spec successively broadened by different velocity dispersions in [km/s]. The code is not optimized for cpu speed, but it shows the principle of how it works.

<img src="https://github.com/neutrinomuon/PyKinematicalBroadening/blob/main/figures/KinematicalBroadening.png?raw=true" width="90%">

#### <b>LICENSE</b>

Attribution-NonCommercial-NoDerivatives 4.0 (CC BY-NC-ND 4.0)

<img src="https://github.com/neutrinomuon/PyKinematicalBroadening/blob/main/figures/cc_logo.png?raw=true" width="10%">

<a href='https://creativecommons.org/licenses/by-nc-nd/4.0/'>Creative Commons Attribution-NonCommercial-NoDerivs (CC-BY-NC-ND)</a>

