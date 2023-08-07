.. SRT Single Dish Tools documentation master file, created by
   sphinx-quickstart on Tue Jan 19 18:32:56 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the SRT Single Dish Tools documentation!
===================================================

Introduction
------------

The Sardinia Radio Telescope Single Dish Tools (SDT) are a set of Python (>=3.6)
tools designed for the quicklook and analysis of single-dish radio data,
starting from the backends present at the Sardinia Radio Telescope.
They are composed of a Python 3.4+ library for developers
and a set of command-line scripts to soften the learning curve for new users.

The Python library is written following the modern coding standards
documented in the Astropy Coding Guidelines. Automatic tests cover
a significant fraction of the code, and are launched each time a commit
is pushed to the `Github`_ repository, through Github Actions.
The current version is 0.7-devXXX, indicating the development version towards 0.7.
See below the tentative roadmap.

In the current implementation, spectroscopic and total-power on-the-fly
scans are supported, both as part of standalone flux measurements through
"cross scans" and as parts
of a map. Maps are formed through a series of scans that swipe the source
region.

.. figure:: images/otf_vs_xsc.jpg
   :width: 80 %
   :alt: otf vs xscan
   :align: center

   **Figure 1.** On-the-fly maps vs cross scan strategies for single dish
   observations.
   The first is able to produce images, the second is used to obtain quick
   flux measurements of point-like sources.

.. _Github: https://github.com/discos/srt-single-dish-tools

Tentative Roadmap
-----------------

+ v.0.1: Simple map creation, draft calibrated fluxes

+ v.0.2: Stable calibrated fluxes, use of multibeam in the K band

+ v.0.3: Stabilization of interactive interface

+ v.0.4: Generalized, user-supplied scanset filters

+ v.0.5: Converters to MBFITS and CLASS

+ v.0.6: **Sun Images**

+ v.0.7: Full support of general coordinate systems, including Galactic

+ v.0.8: Improved calibration, accept input gain curves

+ v.1.0: Code release.


Installation
------------

Prerequisites
~~~~~~~~~~~~~

Anaconda and virtual environment (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We strongly suggest to install the
`Anaconda <https://www.anaconda.com/download>`__ Python distribution.
Once the installation has finished, you should have a working ``conda``
command in your shell. First of all, create a new environment:

.. code-block:: console

    $ conda create -n py3 python=3

load the new environment:

.. code-block:: console

    $ conda activate py3

and install the dependencies (including a few optional but recommended):

.. code-block:: console

    (py3) $ conda install astropy>=3.6 scipy numpy matplotlib pyyaml h5py statsmodels numba

.. code-block:: console

    $ pip install pyregion


Other Python distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^

Install the dependencies with pip (including a few optional but
recommended):

.. code-block:: console

    $ pip install astropy>=3.6 scipy numpy matplotlib pyyaml h5py statsmodels numba pyregion


Cloning and installation
~~~~~~~~~~~~~~~~~~~~~~~~

Clone the repository:

.. code-block:: console

    (py3) $ cd /my/software/directory/
    (py3) $ git clone https://github.com/discos/srt-single-dish-tools.git

or if you have deployed your SSH key to Github:

.. code-block:: console

    (py3) $ git clone git@github.com:discos/srt-single-dish-tools.git

Then:

.. code-block:: console

    (py3) $ cd srt-single-dish-tools
    (py3) $ pip install .

That's it. After installation has ended, you can verify that software is
installed by executing:

.. code-block:: console

    (py3) $ SDTimage -h

If the help message appears, you're done!

Updating
~~~~~~~~

To update the code, simply run ``git pull`` and reinstall:

.. code-block:: console

    (py3) $ git pull
    (py3) $ pip install .


Tutorials
~~~~~~~~~

.. toctree::
  :maxdepth: 2

  tutorials/imaging
  tutorials/converters

Command line interface
----------------------

.. toctree::
  :maxdepth: 2

  scripts/cli

API documentation
-----------------

.. toctree::
  :maxdepth: 2

  srttools/modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
