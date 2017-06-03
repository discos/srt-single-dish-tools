# SRT Single dish tools #
[![Travis CI status](https://travis-ci.org/matteobachetti/srt-single-dish-tools.svg?branch=master)](https://travis-ci.org/matteobachetti/srt-single-dish-tools)
[![Appveyor build status](https://ci.appveyor.com/api/projects/status/m7ll5da3w8ukkrjd/branch/master?svg=true)](https://ci.appveyor.com/project/matteobachetti/srt-single-dish-tools/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/matteobachetti/srt-single-dish-tools/badge.svg?branch=master)](https://coveralls.io/github/matteobachetti/srt-single-dish-tools?branch=master)
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

    $ source activate py3

and install the dependencies (including a few optional but recommended):

    (py3) $ conda install astropy scipy numpy matplotlib pyyaml h5py statsmodels numba

#### Other Python distributions
Install the dependencies with pip (including a few optional but recommended):

    $ pip install astropy scipy numpy matplotlib pyyaml h5py statsmodels numba

### Cloning and installation

Clone the repository:

    (py3) $ cd /my/software/directory/
    (py3) $ git clone https://github.com/matteobachetti/srt-single-dish-tools.git

or if you have deployed your SSH key to Gitlab:

    (py3) $ git clone git@github.com:matteobachetti/srt-single-dish-tools.git

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

### Contribution guidelines ###

[Why writing contribution guidelines, when Astropy has made such a good job already?](http://docs.astropy.org/en/stable/development/codeguide.html)
