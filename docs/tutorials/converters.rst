Format converters
-----------------
Italian antennas baded on the ACS control system save raw data in a FITS format called ``fitszilla``.
Users of other facilities might find it useful to have data converted in a known format.
The SRT single dish tools have a convenient script for that, called ``SDTconvert``.

.. note::

    If you are using Anaconda: remember to initialize the environment. E.g. for an environment called ``py3``

    .. code-block:: console

        $ source activate py3


Read logs and update information in fits files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Some information is not available in fits files produced with old versions of Nuraghe.
Provided that an ``acs.xml`` file is available in the system, we can get a table of the full log in human-readable format with

.. code-block:: console

    (py3) $ SDTparselog acs.xml --to-csv

There are a few shortcuts to the information we want to retrieve.
E.g. we can get a list of FITS files with the calibration mark on:

.. code-block:: console

    (py3) $ SDTparselog acs.xml --list-calon

The calibration mark ON is one of the missing pieces of information in some versions of Nuraghe.

We can get the list of files where the calibration mark was on like above, and then apply it with ``SDTbulkchange``, with the following shortcut:

.. code-block:: console

    (py3) $ SDTparselog acs.xml --list-calon | grep fits | xargs SDTbulkchange --apply-cal-mark


CLASS format
~~~~~~~~~~~~
To get the data in a calibrated CLASS format readable into GILDAS, provided that the observations had a compatible ON-OFF or ON-OFF-CAL sequence (for CAL, apply the information to fits files as described above), type

.. code-block:: console

    (py3) $ SDTconvert -f classfits path/to/directory_of_observation

This will save the calibrated data into a directory called ``path/to/directory_of_observation_classfits/``.
To save in the current directory regardless of the path to the original data, use the option ``--save-locally``.
.. code-block:: console

    (py3) $ SDTconvert -f classfits path/to/directory_of_observation --save-locally

This will save the calibrated data into ``./directory_of_observation_classfits/``

We use the FITS format readable into CLASS, and for convenience we also save a small script that, launched from the user's version of CLASS, is able to convert the data into the native CLASS format.
We do not make the direct conversion to the binary CLASS format for portability issues, but we found that in practice the FITS format is understood correctly across the last four years of GILDAS versions.

Simple feed coordinate conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The FITS format used at the SRT only saves the coordinates of the central feed, and the coordinates of the other feeds need to be calculated based on their offsets in the focal plane.

SDT knows how to treat this problem. However, users wanting to analyze the data with their own software can use ``SDTconvert``:

.. code-block:: console

    (py3) $ SDTconvert -f fitsmod directory_of_observation

This will create a separate extension called ``COORD`` *n* for each feed, where *n* is the number of the feed. Feed 0 will not need a separate extension. Each extension will contain the updated right ascension and declination of the sky region observed by each feed.

MBFITS
~~~~~~
Many European facilities use MBFITS as their raw data format.
``SDTconvert`` can convert the raw data from Italian facilities to this format.

To get the data in the Hierarchical MBFITS format, with the scan divided in multiple files under a directory tree, use

.. code-block:: console

    (py3) $ SDTconvert -f mbfits directory_of_observation

To get a single MBFITS file for each Frontend-Backend combination, use instead

.. code-block:: console

    (py3) $ SDTconvert -f mbfitsw directory_of_observation

SDFITS
~~~~~~
CASA and other software accept data in the SDFITS format.

To get the data in the SDFITS format, with the scan divided in multiple files under a directory tree, use

.. code-block:: console

    (py3) $ SDTconvert -f sdfits directory_of_observation

