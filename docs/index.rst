.. SRT Single Dish Tools documentation master file, created by
   sphinx-quickstart on Tue Jan 19 18:32:56 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SRT Single Dish Tools's documentation!
=================================================

Installation
------------

Anaconda and virtual environment (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We strongly suggest to install the
`Anaconda <https://www.continuum.io/downloads>`__ Python distribution.
Once the installation has finished, you should have a working ``conda``
command in your shell. First of all, create a new environment:

::

    $ conda create -n py3 python=3

load the new environment:

::

    $ source activate py3

and install the dependencies (including a few optional but recommended):

::

    (py3) $ conda install astropy scipy numpy matplotlib pyyaml h5py statsmodels numba


Other Python distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^

Install the dependencies with pip (including a few optional but
recommended):

::

    $ pip install astropy scipy numpy matplotlib pyyaml h5py statsmodels numba

Cloning and installation
~~~~~~~~~~~~~~~~~~~~~~~~

Clone the repository:

::

    (py3) $ cd /my/software/directory/
    (py3) $ git clone https://github.com/matteobachetti/srt-single-dish-tools.git

or if you have deployed your SSH key to Github:

::

    (py3) $ git clone git@github.com:matteobachetti/srt-single-dish-tools.git

Then:

::

    (py3) $ cd srt-single-dish-tools
    (py3) $ python setup.py install

That's it. After installation has ended, you can verify that software is
installed by executing:

::

    (py3) $ SDTimage -h

If the help message appears, you're done!

Updating
~~~~~~~~

To update the code, simply run ``git pull`` and reinstall:

::

    (py3) $ git pull
    (py3) $ python setup.py install


Quick introduction
------------------

Calibrated light curves
^^^^^^^^^^^^^^^^^^^^^^^
Go to a directory close to your data set. For example::

    (py3) $ ls
    observation1/
    observation2/
    calibrator1/
    calibrator2/
    observation3/
    calibrator3/
    calibrator4/
    (....)

It is not required that scan files are directly inside ``observation1`` etc.,
they might be inside subdirectories. The important thing is to correctly point
to them in the configuration file as explained below.

Produce a dummy calibration file, to be modified, with::

    (py3) $ SDTlcurve --sample-config

This produces a boilerplate configuration file, that we modify to point to our
observations, and give the correct information to our program::

    (py3) $ mv sample_config_file.ini MySource.ini  # give a meaningful name!
    (py3) $ emacs MySource.ini

    (... modify file...)

    (py3) $ cat sample_config_file.ini
    (...)
    [analysis]
    (...)
    list_of_directories :
    ;;Two options: either a list of directories:
        dir1
        dir2
        dir3
    calibrator_directories :
        cal1
        cal2
    noise_threshold : 5

    ;; Channels to save from RFI filtering. It might indicate known strong spectral
    ;; lines
    goodchans :

Finally, execute the light curve creation. If data were taken with a Total
Power-like instrument and they do not contain spectral information, it is
sufficient to run::

    (py3) $ SDTlcurve -c MySource.ini

Otherwise, specify the minimum and maximum frequency to select in the spectrum,
with the ``--splat`` option::

    (py3) $ SDTlcurve -c MySource.ini --splat <freqmin>:<freqmax>

where ``freqmin``, ``freqmax`` are in MHz referred to the *minimum* frequency
of the interval. E.g. if our local oscillator is at 6900 MHz and we want to cut
from 7000 to 7500, ``freqmin`` and ``freqmax`` will be 100 and 600 resp.
The above command will:

+ Run through all the scans in the directories specified in the config file

+ Clean them up with a rough but functional algorithm for RFI removal that makes use of the spectral information

+ Create a csv file for each source, containing three columns: time, flux, flux error for each cross scan

The light curve will also be saved in a text file.

Images from OTF maps
~~~~~~~~~~~~~~~~~~~~
The procedure is mostly the same as for light curves.
Go to a directory close to your data set. For example::

    (py3) $ ls
    observation1/
    observation2/
    observation3/
    (....)

It is not required that scan files are directly inside ``observation1`` etc.,
they might be inside subdirectories. The important thing is to correctly point
to them in the configuration file as explained below.

Produce a dummy calibration file, to be modified, with::

    (py3) $ SDTimage --sample-config

This produces a boilerplate configuration file. Give it a meaningful name::

    (py3) $ mv sample_config_file.ini MySource.ini

Modify the file point to our observations, and give the correct information to
our program. Follow the same indentation as in the examples. Comments are done
with a semicolon.
Pay particular attention to the `pixel_size` and to the `list_of_directories`,
pointing to the directories containing scans.::

    (py3) $ emacs MySource.ini

    (... modify file...)

    (py3) $ cat MySource.ini
    [local]
    ; the directory where the analysis will be executed.
        workdir : .
    ; the root directory of the data repository.
        datadir : .

    [analysis]
        projection : ARC
        interpolation : spline
        prefix : test_
        list_of_directories :
            observation1/
            observation2/
            morestuff/observation3/

    ;; Coordinates have to be specified in decimal degrees. ONLY use if different
    ;; from target coordinates!
    ;    reference_ra : 10.5
    ;    reference_dec : 5.3

    ;; Pixel size in arcminutes

        pixel_size : 1

    ;; Channels to save from RFI filtering. It might indicate known strong spectral
    ;; lines
        goodchans :

        noise_threshold : 10


Finally, execute the map calculation. If data were taken with a Total
Power-like instrument and they do not contain spectral information, it is
sufficient to run::

    (py3) $ SDTimage -c MySource.ini

Otherwise, specify the minimum and maximum frequency to select in the spectrum,
with the ``--splat`` option::

    (py3) $ SDTimage -c MySource.ini --splat <freqmin>:<freqmax>

where ``freqmin``, ``freqmax`` are in MHz referred to the *minimum* frequency
of the interval. E.g. if our local oscillator is at 6900 MHz and we want to cut
from 7000 to 7500, ``freqmin`` and ``freqmax`` will be 100 and 600 resp.
The above command will:

+ Run through all the scans in the directories specified in the config file

+ Clean them up with a rough but functional algorithm for RFI removal that makes use of the spectral information

+ Create a single frequency channel per polarization by summing the contributions between ``freqmin`` and ``freqmax``, and discarding the rest;

+ Create the map in FITS format readable by DS9. The FITS extensions IMGCH0, IMGCH1, etc. contain an image for each polarization channel.

Advanced imaging (TBC)
~~~~~~~~~~~~~~~~~~~~~~
The automatic RFI removal procedure is often unable to clean all the data.
The map might have some residual "stripes" due to bad scans. No worries! Launch
the above command with the ``--interactive`` option::

    (py3) $ SDTimage -c MySource.ini --splat <freqmin>:<freqmax> --interactive

This will open a screen like this:

    <placeholder>

where on the right you have the current status of the image, and on the left,
larger, an image of the *standard deviation* of the pixels. Pixels with higher
standard deviation might be due to a real source with high variability or high
flux gradients, or to interferences. On this standard deviation image, you can
point with the mouse and press 'A' on the keyboard to load all scans passing
through that pixel. A second window will appear with a bunch of scans.

    <placeholder>

Click on a bad scan and filter it according to the instructions printed in the
terminal.

Calibration of images
~~~~~~~~~~~~~~~~~~~~~
First of all, call::

    (py3) $ SDTcal  --sample-config

Modify the configuration file adding calibrator directories below `calibrator_directories`::

   calibrator_directories :
      datestring1-3C295/
      datestring2-3C295/

Then, call again ``SDTcal`` with the ``--splat`` option, using **the same frequency range**
of the sources.::

    (py3) $ SDTcal -c MyCalibrators.ini --splat <freqmin>:<freqmax> -o calibration.hdf5

Then, call ``SDTimage`` with the ``--calibrate`` option, as follows::

    (py3) $ SDTimage --calibrate calibration.hdf5 -c MySource.ini --splat <freqmin>:<freqmax> --interactive

... and that's it!



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
