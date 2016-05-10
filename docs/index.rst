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
    (py35) $ git clone https://username@bitbucket.org/mbachett/srt-single-dish-tools.git

or if you have deployed your SSH key to Bitbucket::

    (py35) $ git clone git@bitbucket.org:mbachett/srt-single-dish-tools.git

Then::

    (py35) $ cd srt-single-dish-tools
    (py35) $ python setup.py install

That's it. After installation has ended, you can verify that software is installed
by executing::

    (py35) $ SDTlcurve -h

If the help message appears, you're done!

Updating
~~~~~~~~

To update the code, simply run ``git pull`` and reinstall::

    (py35) $ git pull
    (py35) $ python setup.py install

Quick introduction
------------------
Calibrated light curves (OUTDATED)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Go to a directory close to your data set. For example::

    (py35) $ ls
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

    (py35) $ SDTlcurve --sample-config

This produces a boilerplate configuration file, that we modify to point to our
observations, and give the correct information to our program::

    (py35) $ mv sample_config_file.ini MySource.ini  # give a meaningful name!
    (py35) $ emacs MySource.ini

    (... modify file...)

    (py35) $ cat sample_config_file.ini
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
            calibrator1/
            calibrator2/
            observation3/
            calibrator3/
            calibrator4/
    ;; Channels to save from RFI filtering. It might indicate known strong spectral
    ;; lines
        goodchans :

Finally, execute the light curve creation. If data were taken with a Total
Power-like instrument and they do not contain spectral information, it is
sufficient to run::

    (py35) $ SDTlcurve -c MySource.ini

Otherwise, specify the minimum and maximum frequency to select in the spectrum,
with the ``--splat`` option::

    (py35) $ SDTlcurve -c MySource.ini --splat <freqmin>:<freqmax>

where ``freqmin``, ``freqmax`` are in MHz referred to the *minimum* frequency
of the interval. E.g. if our local oscillator is at 6900 MHz and we want to cut
from 7000 to 7500, ``freqmin`` and ``freqmax`` will be 100 and 600 resp.
The above command will:

+ Run through all the scans in the directories specified in the config file

+ Clean them up with a rough but functional algorithm for RFI removal that makes use of the spectral information

+ Diplay the following products:

  1.  The calibrated light curve with statistical error bars arising from the fit + a grey band indicating the systematic error that might arise from the normalization of the tabulated calibrator fluxes.

  2. Plots of counts-to-Jansky conversion versus elevation: there will often be a linear trend here.

  3. Plots of source misalignment w.r.t. elevation: this will often be surprisingly high

  4. All fitted scans, on a per-source basis.

The light curve will also be saved in a text file.

Images from OTF maps
~~~~~~~~~~~~~~~~~~~~
The procedure is mostly the same as for light curves.
Go to a directory close to your data set. For example::

    (py35) $ ls
    observation1/
    observation2/
    observation3/
    (....)

It is not required that scan files are directly inside ``observation1`` etc.,
they might be inside subdirectories. The important thing is to correctly point
to them in the configuration file as explained below.

Produce a dummy calibration file, to be modified, with::

    (py35) $ SDTimage --sample-config

This produces a boilerplate configuration file. Give it a meaningful name::

    (py35) $ mv sample_config_file.ini MySource.ini

Modify the file point to our observations, and give the correct information to
our program. Follow the same indentation as in the examples. Comments are done
with a semicolon.
Pay particular attention to the `pixel_size` and to the `list_of_directories`,
pointing to the directories containing scans.::

    (py35) $ emacs MySource.ini

    (... modify file...)

    (py35) $ cat MySource.ini
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

    (py35) $ SDTimage -c MySource.ini

Otherwise, specify the minimum and maximum frequency to select in the spectrum,
with the ``--splat`` option::

    (py35) $ SDTimage -c MySource.ini --splat <freqmin>:<freqmax>

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

    (py35) $ SDTimage -c MySource.ini --splat <freqmin>:<freqmax> --interactive

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

    (py35) $ SDTcal  --sample-config

Modify the configuration file adding calibrator directories below `calibrator_directories`::

   calibrator_directories :
      datestring1-3C295/
      datestring2-3C295/

Then, call again ``SDTcal`` with the ``--splat`` option, using **the same frequency range**
of the sources.::

    (py35) $ SDTcal -c MyCalibrators.ini --splat <freqmin>:<freqmax> -o calibration.hdf5

Then, call ``SDTimage`` with the ``--calibrate`` option, as follows::

    (py35) $ SDTimage --calibrate calibration.hdf5 -c MySource.ini --splat <freqmin>:<freqmax> --interactive

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
