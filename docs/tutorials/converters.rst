Format converters
-----------------
Italian antennas baded on the ACS control system save raw data in a FITS format called ``fitszilla``.
Users of other facilities might find it useful to have data converted in a known format.
The SRT single dish tools have a convenient script for that, called ``SDTconvert``.


CLASS format
~~~~~~~~~~~~
To get the data in a calibrated CLASS format readable into GILDAS, provided that the observations had a compatible ON-OFF or ON-OFF-CAL sequence, type

.. code-block:: console

    (py3) $ SDTconvert -f classfits directory_of_observation

This will save the calibrated data into a directory.
We use the FITS format readable into CLASS, and for convenience we also save a small script that, launched from the user's version of CLASS, is able to convert the data into the native CLASS format.
We do not make the direct conversion to the binary CLASS format for portability issues, but we found that in practice the FITS format is understood correctly across the last four years of GILDAS versions.

Simple feed coordinate conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The FITS format used at the SRT only saves the coordinates of the central feed, and the coordinates of the other feeds need to be calculated based on their offsets in the focal plane.

SDT knows how to treat this problem. However, users wanting to analyze the data with their own software can use ``SDTconvert``:

.. code-block:: console

    (py3) $ SDTconvert -f fitsmod directory_of_observation

This will create a separate extension called ``COORD``*n* for each feed, where *n* is the number of the feed. Feed 0 will not need a separate extension. Each extension will contain the updated right ascension and declination of the sky region observed by each feed.

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

