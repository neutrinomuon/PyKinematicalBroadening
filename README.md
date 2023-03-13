### PyKinematicalBroadening
email: [antineutrinomuon@gmail.com](mailto:antineutrinomuon@gmail.com), [jean@astro.up.pt](mailto:jean@astro.up.pt)

© Copyright ®

J.G. - Jean Gomes

<hr>

![python3](https://img.shields.io/pypi/pyversions/pyfluxconserving)

<hr>

#### <b>RESUME</b>

<img src="https://raw.githubusercontent.com/neutrinomuon/PyKinematicalBroadening/main/figures/Kinematics.png" width=120>

PyKinematicalBroadening is a Python repository for Extragalactic Kinematics that performs kinematical broadening of a spectrum by convolving it with a kernel in velocity space, which results in a broadened model. This code provides a function called 'broadening,' which uses a Gaussian kernel for convolution. However, in addition to the Gaussian kernel, the code can also use Gauss-Hermite kernels for convolution.

The Gauss-Hermite kernels provide a more general perspective than the Gaussian kernels, as they can account for higher-order velocity moments of the broadening process. The code defines the width and mean velocity of the kernel with vd_sigma and vc0_gals, respectively, while the number of points for the kernel is set with the Ni_Gauss parameter. The output spectrum is defined at wavelengths lambda_s and is returned as fluxes_s. The fill_val parameter defines the value to use for regions outside of the original wavelength range, and verbosity controls the level of detail of console output.

In summary, the PyKinematicalBroadening code provides a flexible way to apply kinematical broadening to a spectrum using Gaussian or Gauss-Hermite kernels. The code reads in a test spectrum, interpolates it onto a set of equally spaced wavelength values, and loops over different velocity dispersions to call the 'broadening' function, which broadens the spectrum and plots the results.

<hr>

#### <b>INSTALLATION</b>

You can easily install <a href=https://pypi.org/project/PyKinematicalBroadening/>PyKinematicalBroadening</a> by using pip - <a href='https://pypi.org/'>PyPI - The Python Package Index</a>:

<pre>
pip install PyKinematicalBroadening
</pre>

<br>or by using a generated conda repository <a href='https://anaconda.org/neutrinomuon/PyKinematicalBroadening'>https://anaconda.org/neutrinomuon/PyKinematicalBroadening</a>:

[![badgetanaconda](https://anaconda.org/neutrinomuon/PyKinematicalBroadening/badges/version.svg)](https://anaconda.org/neutrinomuon/PyKinematicalBroadening/badges/version.svg)
[![badgetreleasedate](https://anaconda.org/neutrinomuon/PyKinematicalBroadening/badges/latest_release_date.svg)](https://anaconda.org/neutrinomuon/PyKinematicalBroadening/badges/latest_release_date.svg)
[![badgetplatforms](https://anaconda.org/neutrinomuon/PyKinematicalBroadening/badges/platforms.svg
)](https://anaconda.org/neutrinomuon/PyKinematicalBroadening/badges/platforms.svg)

<pre>
conda install -c neutrinomuon PyKinematicalBroadening
</pre>

<br>OBS.: Linux, OS-X ad Windows pre-compilations available in conda.

You can also clone the repository and install by yourself in your machine:

<pre>
git clone https://github.com/neutrinomuon/PyKinematicalBroadening
python setup.py install
</pre>

<hr>

#### <b>EXAMPLE</b>

Example of the test_spectrum test_spectrum.spec successively broadened by different velocity dispersions in [km/s]. The code is not optimized for cpu speed, but it shows the principle of how it works.

<img src="https://github.com/neutrinomuon/PyKinematicalBroadening/blob/main/figures/KinematicalBroadening.png?raw=true" width="90%">

<hr>

#### <b>STRUCTURE</b>

<pre>
PyKinematicalBroadening
├── MANIFEST.in
├── dist
│   ├── PyKinematicalBroadening-0.0.3.tar.gz
│   ├── PyKinematicalBroadening-0.0.5.tar.gz
│   ├── PyKinematicalBroadening-0.0.6.tar.gz
│   └── PyKinematicalBroadening-0.0.4.tar.gz
├── README.md
├── figures
│   ├── KinematicalBroadening.png
│   └── cc_logo.png
├── PyKinematicalBroadening.egg-info
│   ├── PKG-INFO
│   ├── dependency_links.txt
│   ├── SOURCES.txt
│   ├── top_level.txt
│   └── requires.txt
├── LICENSE.txt
├── setup.py
├── tutorials
│   ├── .ipynb_checkpoints
│   │   └── Example 1 - Kinematical Broadening-checkpoint.ipynb
│   └── Example 1 - Kinematical Broadening.ipynb
├── pykinematicalbroadening
│   ├── win-32
│   │   └── pykinematicalbroadening-0.0.5-py39_0.tar.bz2
│   ├── linux-armv7l
│   │   └── pykinematicalbroadening-0.0.5-py39_0.tar.bz2
│   ├── linux-armv6l
│   │   ├── .projectignore
│   │   └── pykinematicalbroadening-0.0.5-py39_0.tar.bz2
│   ├── linux-s390x
│   │   └── pykinematicalbroadening-0.0.5-py39_0.tar.bz2
│   ├── linux-ppc64
│   │   └── pykinematicalbroadening-0.0.5-py39_0.tar.bz2
│   ├── linux-aarch64
│   │   ├── .projectignore
│   │   └── pykinematicalbroadening-0.0.5-py39_0.tar.bz2
│   ├── linux-32
│   │   ├── .projectignore
│   │   └── pykinematicalbroadening-0.0.5-py39_0.tar.bz2
│   ├── linux-64
│   │   ├── .projectignore
│   │   └── pykinematicalbroadening-0.0.5-py39_0.tar.bz2
│   ├── osx-64
│   │   └── pykinematicalbroadening-0.0.5-py39_0.tar.bz2
│   ├── meta.yaml
│   ├── win-64
│   │   └── pykinematicalbroadening-0.0.5-py39_0.tar.bz2
│   ├── README.txt
│   ├── linux-ppc64le
│   │   └── pykinematicalbroadening-0.0.5-py39_0.tar.bz2
│   └── osx-arm64
│       └── pykinematicalbroadening-0.0.5-py39_0.tar.bz2
├── Pykinematicalbroadening.egg-info
│   ├── PKG-INFO
│   ├── dependency_links.txt
│   ├── SOURCES.txt
│   ├── top_level.txt
│   └── requires.txt
├── src
│   └── python
│       ├── __pycache__
│       ├── test_spectrum.spec
│       ├── __init__.py
│       └── PyKinematicalBroadening.py
├── version.txt
└── build
    └── lib
        ├── Pykinematicalbroadening
        └── PyKinematicalBroadening

26 directories, 44 files
</pre>

<hr>

#### <b>REFERENCES</b>

<ol>

<il>Tonry, J., & Davis, M. 1979, AJ, 84, 1511</il>

<il>van der Marel, R. P., & Franx, M. 1993, ApJ, 407, 525</il>
        
<il>Emsellem, E., et al. "The SAURON project - III. Integral-field absorption-line
kinematics of 48 elliptical and lenticular galaxies." Monthly Notices of the Royal Astronomical Society, Volume 352, Issue 3, pp. 721-743. DOI:  
<a href="https://doi.org/10.1111/j.1365-2966.2004.07948.x">10.1111/j.1365-2966.2004.07948.x</a> ; <a href="https://doi.org/10.48550/arXiv.astro-ph/0404034">10.48550/arXiv.astro-ph/0404034</a>. Available
at: <a href="https://ui.adsabs.harvard.edu/abs/2004MNRAS.352..721E/abstract">https://ui.adsabs.harvard.edu/abs/2004MNRAS.352..721E/abstract</a>.</il>

[//]: # (<il>Faber, S. M. "The Stellar Population Histories of Elliptical Galaxies: A
[//]: # Review." Annual Review of Astronomy and Astrophysics, vol. 46, no. 1, 2008,
[//]: # pp. 121-157. DOI: <a
[//]: # href="https://doi.org/10.1146/annurev-astro-082708-101650">10.1146/annurev-astro-082708-101650</a>. Available
[//]: # at: <a
[//]: # href="https://www.annualreviews.org/doi/10.1146/annurev-astro-082708-101650">https://www.annualreviews.org/doi/10.1146/annurev-astro-082708-
[//]: # 101650</a>.</il>)

[//]: # (<il>Peletier, R. F., et al. "The SAURON project - XI. Stellar populations from
[//]: # absorption-line strength maps of 24 early-type spirals." Monthly Notices of
[//]: # the Royal Astronomical Society, vol. 379, no. 2, 2007, pp. 445-469. DOI: <a
[//]: # href="https://doi.org/10.1111/j.1365-2966.2007.11803.x">10.1111/j.1365-2966.2007.11803.x</a>. Available
[//]: # at: <a
[//]: # href="https://academic.oup.com/mnras/article/379/2/445/1078958">https://academic.oup.com/mnras/article/379/2/445/1078958</a>.</il>)

[//]: # (<il>Maraston, C. "Spectral Synthesis of Stellar Populations with Star Formation
[//]: # Histories." Monthly Notices of the Royal Astronomical Society, vol. 362,
[//]: # no. 3, 2005, pp. 799-825. DOI: <a
[//]: # href="https://doi.org/10.1111/j.1365-2966.2005.09340.x">10.1111/j.1365-2966.2005.09340.x</a>. Available
[//]: # at: <a
[//]: # href="https://academic.oup.com/mnras/article/362/3/799/986891">https://academic.oup.com/mnras/article/362/3/799/986891</a>.</il>)

</ol>

<hr>

#### <b>LICENSE</b>

Attribution-NonCommercial-NoDerivatives 4.0 (CC BY-NC-ND 4.0)

<img src="https://github.com/neutrinomuon/PyKinematicalBroadening/blob/main/figures/cc_logo.png?raw=true" width="10%">

<a href='https://creativecommons.org/licenses/by-nc-nd/4.0/'>Creative Commons Attribution-NonCommercial-NoDerivs (CC-BY-NC-ND)</a>

