# SRT Single dish tools #
[![CI Tests](https://github.com/discos/srt-single-dish-tools/workflows/CI%20Tests/badge.svg)](https://github.com/discos/srt-single-dish-tools)
[![codecov](https://img.shields.io/codecov/c/github/discos/srt-single-dish-tools/master.svg?maxAge=0)](https://codecov.io/gh/discos/srt-single-dish-tools)
[![Powered by Astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
[![Documentation Status](https://readthedocs.org/projects/srt-single-dish-tools/badge/?version=latest)](http://srt-single-dish-tools.readthedocs.io/en/latest/?badge=latest)

## Installation

### Preparation and dependencies

#### Anaconda and virtual environment (recommended)
We strongly suggest to install the
[Anaconda](https://www.continuum.io/downloads) Python distribution.
Once the installation has finished, you should have a working `conda`
command in your shell. First of all, create a new environment:

    $ conda create -n py3 python=3

load the new environment:

    $ conda activate py3

and install the dependencies (including a few optional but recommended):

    (py3) $ conda install astropy scipy numpy matplotlib pyyaml h5py statsmodels numba

#### Other Python distributions
Install the dependencies with pip (including a few optional but recommended):

    $ pip install astropy scipy numpy matplotlib pyyaml h5py statsmodels numba

### Cloning and installation

Clone the repository:

    (py3) $ cd /my/software/directory/
    (py3) $ git clone https://github.com/discos/srt-single-dish-tools.git

or if you have deployed your SSH key to Github:

    (py3) $ git clone git@github.com:discos/srt-single-dish-tools.git

Then:

    (py3) $ cd srt-single-dish-tools
    (py3) $ python setup.py install

That's it. After installation has ended, you can verify that software is
installed by executing:

    (py3) $ SDTimage -h

If the help message appears, you're done!

### Updating

To update the code, simply run `git pull` and reinstall:

    (py3) $ git pull
    (py3) $ python setup.py install

### Contribution guidelines

Please follow the [Astropy contribution guidelines](http://docs.astropy.org/en/stable/development/workflow/development_workflow.html), and the [Astropy coding guidelines](http://docs.astropy.org/en/stable/development/codeguide.html#coding-style-conventions). This code is written in Python 3.8+. Tests run at each commit during Pull Requests, so it is easy to single out points in the code that break this compatibility.
See the file Contributing.md for more details

### If you use this code

First of all... **This code is under development!**... so, it might well be that something does not work as expected. For any inquiries, bug reports, or suggestions, please use the [Issues](https://github.com/discos/srt-single-dish-tools/issues) page.

If you used this software package to reduce data for a publication, please write in the acknowledgements something along these lines:

    This work makes use of the [SRT Single Dish Tools](https://github.com/discos/srt-single-dish-tools/) .

We will submit it to the [Astrophysics Source Code Library](www.ascl.net) soon, and update the text accordingly.
