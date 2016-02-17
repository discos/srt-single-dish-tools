.. SRT Single Dish Tools documentation master file, created by
   sphinx-quickstart on Tue Jan 19 18:32:56 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SRT Single Dish Tools's documentation!
=================================================

Installation
------------

Anaconda and virtual environment (recommended but optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We strongly suggest to install the [Anaconda](https://www.continuum.io/downloads)
Python distribution. Once the installation has finished, you should have a working
``conda`` command in your shell. First of all, create a new environment::

    $ conda create -n py35 python=3.5

load the new environment::

    $ source activate py35

and install the dependencies::

    (py35) $ conda install matplotlib h5py astropy scipy numpy

Cloning and installation
~~~~~~~~~~~~~~~~~~~~~~~~

Clone the repository::

    (py35) $ cd /my/software/directory/
    (py35) $ git clone https://username@bitbucket.org/srt-tools/srt-single-dish-tools.git

or if you have deployed your SSH key to Bitbucket::

    (py35) $ git clone git@bitbucket.org:srt-tools/srt-single-dish-tools.git

Then::

    (py35) $ cd srt-single-dish-tools
    (py35) $ python setup.py install

That's it. After installation has ended, you can verify that software is installed
by executing::

    (py35) $ SDTlcurve -h

If the help message appears, you're done!

Command line interface
----------------------

.. toctree::
  :maxdepth: 2

  scripts/cli

API documentation
-----------------

.. toctree::
  :maxdepth: 2

  srttools/core/modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
