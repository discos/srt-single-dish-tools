.. SRT Single Dish Tools documentation master file, created by
   sphinx-quickstart on Tue Jan 19 18:32:56 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the SRT Single Dish Tools documentation!
===================================================

Introduction
------------
The Sardinia Radio Telescope Single Dish Tools (SDT) are a set of Python
tools designed for the quicklook and analysis of single-dish radio data,
starting from the backends present at the Sardinia Radio Telescope.
They are composed of a Python (2.7, 3.4+) library for developers
and a set of command-line scripts to soften the learning curve for new users.

The Python library is written following the modern coding standards
documented in the Astropy Coding Guidelines. Automatic tests cover
a significant fraction of the code, and are launched each time a commit
is pushed to the `Github`_ repository.
The Continuous Integration service `Travis CI`_ is used for that.
The current version is 0.5-devXXX, indicating the development version towards 0.5.
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

.. _Travis CI: http://www.travis-ci.com
.. _Github: https://github.com/matteobachetti/srt-single-dish-tools

Tentative Roadmap
-----------------

+ v.0.1: Simple map creation, draft calibrated fluxes

+ v.0.2: Stable calibrated fluxes, use of multibeam in the K band

+ v.0.3: Stabilization of interactive interface

+ v.0.4: Generalized, user-supplied scanset filters

+ v.**0.5**: Converters to MBFITS and CLASS

+ v.0.6: Improved calibration, accept input gain curves

+ v.0.7: Improved RFI support, using simple techniques of machine learning

+ v.0.8: Full support of general coordinate systems, including Galactic

+ v.1.0: code release.


Installation
------------

Prerequisites
~~~~~~~~~~~~~

Anaconda and virtual environment (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We strongly suggest to install the
`Anaconda <https://www.continuum.io/downloads>`__ Python distribution.
Once the installation has finished, you should have a working ``conda``
command in your shell. First of all, create a new environment:

.. code-block:: console

    $ conda create -n py3 python=3

load the new environment:

.. code-block:: console

    $ source activate py3

and install the dependencies (including a few optional but recommended):

.. code-block:: console

    (py3) $ conda install astropy>=3 scipy numpy matplotlib pyyaml h5py statsmodels numba

.. code-block:: console

    $ pip install pyregion


Other Python distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^

Install the dependencies with pip (including a few optional but
recommended):

.. code-block:: console

    $ pip install astropy>=3 scipy numpy matplotlib pyyaml h5py statsmodels numba pyregion

Cloning and installation
~~~~~~~~~~~~~~~~~~~~~~~~

Clone the repository:

.. code-block:: console

    (py3) $ cd /my/software/directory/
    (py3) $ git clone https://github.com/matteobachetti/srt-single-dish-tools.git

or if you have deployed your SSH key to Github:

.. code-block:: console

    (py3) $ git clone git@github.com:matteobachetti/srt-single-dish-tools.git

Then:

.. code-block:: console

    (py3) $ cd srt-single-dish-tools
    (py3) $ python setup.py install

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
    (py3) $ python setup.py install



Tutorial
--------
In this tutorial, we will see how to obtain calibrated images and
light curves from a set of on-the-fly (OTF) scans done with the SRT.
Data are taken with the SARDARA ROACH2-based backend, with a
bandwidth of 1024 MHz and 1024 channels.
In this tutorial we will first learn how the software does a
semi-automatic cleaning of the data from radio-frequency
interferences (RFI), and how to tweak the relevant parameters to do
the cleaning properly.
Then, we will generate rough images with the default baseline
subtraction algorithms.
Afterwards, we will load a set of calibrators to perform the
conversion from signal level to Janskys/pixel. Finally, we will
apply the calibration to the previously generated images.

Inspect the observation
~~~~~~~~~~~~~~~~~~~~~~~
During a night of observations, we will in general observe a number of calibrators
and sources, in random order. Our observation will be split into a series
of directories:

.. code-block:: console

    (py3) $ ls
    2016-05-04-220022_Src1/
    2016-05-04-223001_Src1/
    2016-05-04-230001_Cal1/
    2016-05-04-230200_Cal2/
    2016-05-04-230432_Src1/
    2016-05-04-233523_Src1/
    (....)

Some of these observations might have been done in different bands, or using different
receivers, and you might have lost the list of observations (or the user was not the observer).
The script `SDTinspector` is there to help, dividing the observations in groups based
on observing time, backend, receiver, etc.:

.. code-block:: console

    (py3) $ SDTinspect */
    Group 0, Backend = ROACH2, Receiver = CCB
    ------------------
    Src1, observation 0

    Source observations:
    2016-05-04-220022_Src1/
    2016-05-04-223001_Src1/
    2016-05-04-230432_Src1/

    Calibrator observations:
    2016-05-04-230001_Cal1/
    2016-05-04-230200_Cal2/

    Group 1, Backend = ROACH2, Receiver = KKG
    ---------------
    Src1, observation 1

    Source observations:
    2016-05-04-233523_Src1/
    (.....)

    Calibrator observations:
    (.....)

With the ``-d`` option, the script will also dump automatically
a set of config files ready for the next step in the analysis:

.. code-block:: console

    (py3) $ SDTinspect */
    Group 0, Backend = ROACH2, Receiver = CCB
    (.....)
    (py3) $ ls -alrt
    CCB_ROACH_Src1_Obs0.ini
    KKG_ROACH_Src1_Obs1.ini

Modify config files
~~~~~~~~~~~~~~~~~~~
If you did not pre-generate config files with the procedure above,
you can generate a boilerplate config file with:

.. code-block:: console

    (py3) $ SDTlcurve --sample-config
    (py3) $ ls
    (...)
    sample_config_file.ini

In the following, we will use the config files generated by SDTinspect,
but it is very easy to adapt to the case of a custom-modified boilerplate.

Config files have this overall structure (slight changes might occur, like
equals signs being changed to semicolons):

.. code-block:: console

    (py3) $ cat CCB_ROACH_Src1_Obs0.ini
    [local]
    workdir = .
    datadir = .

    [analysis]
    projection = ARC
    interpolation = spline
    list_of_directories =
        2016-05-04-220022_Src1/
        2016-05-04-223001_Src1/
        2016-05-04-230432_Src1/
    calibrator_directories =
        2016-05-04-230001_Cal1/
        2016-05-04-230200_Cal2/
    noise_threshold = 5
    pixel_size = 1
    goodchans =

You will likely not change the kind of interpolation or the projection
in the plane of the sky (but if instead of ``ARC`` you want something
different, `all projections in this list are supported`_).
``goodchans`` is a list of channels that can be excluded from
automatic filtering (for example, because they might contain an important
spectral line.)

``pixel_size`` is by default 1 arcminute. You might want to change this
depending on the density of scans and the beam size at the observing frequency.
Usually, 1/3 of the beam size is ok for dense OTF scan campaigns, while
a larger value is better for sparse observations.

Also, you might know already that some observations were bad. In this case,
it's sufficient to take them out of the list above.

.. _all projections in this list are supported: http://docs.astropy.org/en/stable/wcs/#supported-projections


Preprocess the files
~~~~~~~~~~~~~~~~~~~~

.. figure:: images/filtered_scan.jpg
   :width: 80 %
   :alt: scan filtering
   :align: center

   **Figure 2.** Output of the automatic filtering procedure for an OTF scan of a calibrator.
   Channels where the root mean square of the signal is too high or too low are
   automatically filtered out. The threshold is encoded in the ``noise_threshold``
   variable in the config file. This is the number of standard deviations from the median
   r.m.s. in a given interval.
   Optionally the user can choose the frequency interval (blue vertical lines).
   In the two right panels, one can see the scan before and after the cleaning.
   In the right-lower panel, the uncleaned scan is reported in grey to help
   the eye.
   The dynamical spectrum before and after the cleaning is
   shown in the two middle panels, and the effect of the cleaning on the scan
   is shown in the two right panels

This step is optional, because it can be merged with image production.
However, for the sake of this tutorial we will proceed in this way
for simplicity.

The easiest way to preprocess an observation is to call ``SDTpreprocess`` on
a config file. The script will load all files, one by one, and do the following
steps:

1. If the backend is spectroscopic, load each scan and filter out all channels whose
   that are more noisy than a given value of rms during the scan, then merge into
   a single channel. As an option (recommended), the user can specify a frequency
   interval that will be merged, otherwise the full frequency interval is taken: for
   this, one can use the option ``--splat <minf:maxf>`` where ``minf``, ``mmaxf``
   are in MHz referred to the *minimum* frequency of the interval. E.g. if our local
   oscillator is at 6900 MHz and we want to cut from 7000 to 7500, ``minf`` and ``mmaxf``
   will be 100 and 600 resp. This process produces plots like the following:

.. code-block:: console

    (py3) $ SDTimage -c CCB_TP_Src1_Obs0.ini --splat 80:1100

2. The single channels that are produced at step 1, or alternatively the single
   channels of a non-spectroscopic backend, will now be processed by a baseline
   subtraction routine. This routine, by default, applies an Asymmetric Least
   Squares Smoothing (`Eilers and Boelens 2005`_) to find the rough alignment
   of the scan, and then improves it by selecting the data that
   are closer to the baseline and making a standard least-square fit.
   This procedure is very fast and aligns the vast majority
   of scans in a fraction of a second. For more complicated scans, an interactive
   interface is also available, albeit with some portability issues that will be
   solved in future versions (use the ``--interactive`` option).
   It is possible to avoid regions with known strong sources. For now, they need
   to be specified by hand, with the ``-e`` option followed by a valid ds9-compatible
   region file containing *circular* regions in the ``fk5`` frame.

3. The results of the first points are saved as ``HDF5`` files in the same directory
   as the original ``fits`` files. This makes it
   much faster to reload the scans for further use. If the user wants to reprocess
   the files from scratch, he/she needs to delete these files first.

.. _Eilers and Boelens 2005: https://zanran_storage.s3.amazonaws.com/www.science.uva.nl/ContentPages/443199618.pdf

Let's produce some images now!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finally, let us execute the map calculation. If data were taken with a Total
Power-like instrument and they do not contain spectral information, it is
sufficient to run

.. code-block:: console

    (py3) $ SDTimage -c CCB_TP_Src1_Obs0.ini

where CCB_TP_Src1_Obs0.ini should be substituted with the wanted config file.
*This is also valid for spectroscopic scans that have already been preprocessed*

.. code-block:: console

    (py3) $ SDTimage -c CCB_ROACH_Src1_Obs0.ini

Otherwise, if preprocessing were not executed before, specify the minimum and
maximum frequency to select in the spectrum,
with the ``--splat`` option (same as before)

.. code-block:: console

    (py3) $ SDTimage -c CCB_ROACH_Src1_Obs0.ini --splat <freqmin>:<freqmax>

The above command will:

+ Run through all the scans in the directories specified in the config file

+ Clean them up if not already done in a previous step, in the same way of ``SDTpreprocess``

+ Create a single frequency channel per polarization by summing the contributions between
  ``freqmin`` and ``freqmax``, and discarding the remaining frequency channels,
  again if not already done in a previous step;

+ Create the map in FITS format readable by DS9. The FITS extensions IMGCH0, IMGCH1,
  etc. contain an image for each polarization channel. The extensions IMGCH<no>-STD
  will contain the *error images* corresponding to IMGH<no>.

The automatic RFI removal procedure might have missed some problematic scan.
The map might have, therefore, some residual "stripes" due to bad scans or wrong
baseline subtraction.

The first thing to do, in these cases, is to go and look at the scans (by going
through the PDF files produced by the calibration process in each subdirectory)
and check that the noise threshold is appropriate for the level of noise found
in scans.
If it is not, as is often the case, and it is sufficient to re-run ``SDTpreprocess``
with the noise threshold changed in the config file to get a better cleaning
of the data.

But ``SDTimage`` has an additional option to align the scans. It's called *global
baseline subtraction*. This procedure makes a *global* fit (option ``-g``) of all scans in an
image, and tries to find the alignment of each scan that minimizes the *total
rms* of the image. This procedure is only valid if the region that is fit is
consistent with having zero average. This is, of course, not valid if the source
is strong. In this case, together with the global fit option, we need to also
specify a set of regions to neglect. This is done in two ways:

+ through a ds9-compatible region file containing *circular* regions in *image* coordinates

+ through the option ``-e`` followed by multiples of three numbers: X, Y and radius,
  in *image* coordinates (SAOimage ds9 or other
  imaging programs can create regions with these coordinates, one just needs to
  copy the numbers.).

In summary, to use the global fitting and discard the region centered at coordinates
x,y=30,33 with radius 10 pixels, run

.. code-block:: console

    (py3) $ SDTimage -g -e 30 33 10 (...additional options)

.. figure:: images/map.png
   :width: 80 %
   :alt: map
   :align: center

   **Figure 3.** Map produced by ``SDTimage``

Advanced imaging (TBC)
~~~~~~~~~~~~~~~~~~~~~~
The automatic RFI removal procedure is often unable to clean all the data.
The map might have some residual "stripes" due to bad scans. No worries! Launch
the above command with the ``--interactive`` option

.. code-block:: console

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
To calibrate the images, one needs to call ``SDTcal`` with the same config files
used for the images if they were produced with ``SDTinspect``. Otherwise, one
can construct an alternative config file with

.. code-block:: console

    (py3) $ SDTcal  --sample-config

and modify the configuration file adding calibrator directories
below `calibrator_directories`

.. code-block:: console

   calibrator_directories :
      datestring1-3C295/
      datestring2-3C295/

Then, call again ``SDTcal`` with the ``--splat`` option, using **the same frequency range**
of the sources.

.. code-block:: console

    (py3) $ SDTcal -c CCB_ROACH_Src1_Obs0.ini--splat <freqmin>:<freqmax> -o calibration.hdf5

Finally, call ``SDTimage`` with the ``--calibrate`` option, e.g.

.. code-block:: console

    (py3) $ SDTimage --calibrate calibration.hdf5 -c CCB_ROACH_Src1_Obs0.ini --splat <freqmin>:<freqmax> --interactive

... and that's it! The image values will be expressed in Jy instead of counts, so that
applying a region with DS9 and calculating the total flux inside the given region will
return the actual total flux contained in the region.

Calibrated light curves
~~~~~~~~~~~~~~~~~~~~~~~
Go to a directory close to your data set. For example

.. code-block:: console

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

Produce a dummy calibration file, to be modified, with

.. code-block:: console

    (py3) $ SDTlcurve --sample-config

This produces a boilerplate configuration file, that we modify to point to our
observations, and give the correct information to our program

.. code-block:: console

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
sufficient to run

.. code-block:: console

    (py3) $ SDTlcurve -c MySource.ini

Otherwise, specify the minimum and maximum frequency to select in the spectrum,
with the ``--splat`` option

.. code-block:: console

    (py3) $ SDTlcurve -c MySource.ini --splat <freqmin>:<freqmax>

where ``freqmin``, ``freqmax`` are in MHz referred to the *minimum* frequency
of the interval. E.g. if our local oscillator is at 6900 MHz and we want to cut
from 7000 to 7500, ``freqmin`` and ``freqmax`` will be 100 and 600 resp.
The above command will:

+ Run through all the scans in the directories specified in the config file

+ Clean them up with a rough but functional algorithm for RFI removal that makes use of the spectral information

+ Create a csv file for each source, containing three columns: time, flux, flux error for each cross scan

The light curve will also be saved in a text file.


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
