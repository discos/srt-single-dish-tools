Cross-scan Tutorial: fluxes and upper limits
--------------------------------------------
Cross scans (a few subscans along RA and Dec through the source position,
repeated many times) are the most efficient way to measure the flux of a
point source of known position, or to put an upper limit on it.
``SDTcross`` does this from the same data and configuration files used by
``SDTimage``.

Preparing the data
~~~~~~~~~~~~~~~~~~
Prepare a configuration file listing the source and calibrator directories,
as described in the :doc:`imaging` tutorial, and run ``SDTcal`` on the
calibrators:

.. code-block:: console

    (py3) $ SDTcal -c MySource.ini --splat <freqmin:freqmax> -o calibration.hdf5

Then, either let ``SDTcross`` read the scans directly from the configuration
file:

.. code-block:: console

    (py3) $ SDTcross -c MySource.ini --calibrate calibration.hdf5

or use a ``_dump.hdf5`` file saved by ``SDTimage`` (e.g. after cleaning the
data interactively):

.. code-block:: console

    (py3) $ SDTimage -c MySource.ini
    (py3) $ SDTcross MySource_dump.hdf5 --calibrate calibration.hdf5

If the dump file already contains calibrated ``<channel>_cal`` columns,
``--calibrate`` can be omitted; without any calibration, fluxes are in
counts.

What SDTcross does
~~~~~~~~~~~~~~~~~~
For each subscan and channel, a smooth baseline (a cubic spline with knots
a few beams apart) and a Gaussian beam centred at the nominal source
position are fitted together. Position and width of the Gaussian are fixed
(a *forced fit*): only its amplitude is measured, and it can be negative
when there is no source.

The flux is the average of the amplitudes over the subscans (LCP and RCP
are averaged subscan by subscan first). Its error comes from the scatter of
the amplitudes between subscans: it is correct also when the noise is not
white (1/f noise, baseline drifts), which is the usual case. Subscans whose
noise is anomalous (e.g. residual RFI) are ignored, with a warning
(``--rfi-nsigma``).

A source is *detected* if its flux is more than ``--nsigma`` (default 3)
times its error. The upper limit is the measured flux (or zero, if negative)
plus ``--nsigma`` times the error, with the calibration error added in
quadrature.

Outputs
~~~~~~~
For an input ``MySource_dump.hdf5`` (or ``MySource.ini``), ``SDTcross``
writes:

- ``MySource_cross.fits``: the main results as keywords of the primary
  header (``FLUX``, ``FLUXERR``, ``UPLIM``, ``NSIGMA``, ``DETECTED``,
  ``BUNIT``, ``CALERR``, the consistency checks ``CONS_ARM``,
  ``CONS_POL`` and ``CONS_SHP``, ...), and three tables: ``SUMMARY`` (results per channel,
  per arm and in total), ``PROFILES`` (the baseline-subtracted profiles) and
  ``SUBSCANS`` (the amplitude in each subscan and channel);
- ``MySource_cross.ecsv``: the summary table;
- ``MySource_cross.png`` and ``MySource_cross_total.png``: the profiles
  along each arm, for each channel and in total, with the fitted beam.

The summary is also printed:

.. code-block:: console

        label      channel  arm   flux      err    nsub    ul    detected
    ------------- --------- --- -------- --------- ---- -------- --------
     Feed0_LCP RA Feed0_LCP  RA 0.008211 0.0007449  523  0.01045     True
    Feed0_LCP Dec Feed0_LCP Dec 0.007619 0.0006957  527 0.009706     True
        Feed0_LCP Feed0_LCP     0.007914 0.0005093 1050 0.009442     True
     Feed0_RCP RA Feed0_RCP  RA 0.006734 0.0005729  518 0.008453     True
    Feed0_RCP Dec Feed0_RCP Dec 0.006997 0.0005473  525 0.008639     True
        Feed0_RCP Feed0_RCP     0.006867 0.0003959 1043 0.008054     True
               RA            RA  0.00754 0.0004811  528 0.008983     True
              Dec           Dec 0.007338 0.0004644  532 0.008731     True
            TOTAL               0.007439 0.0003341 1060 0.008441     True

Checks
~~~~~~
- The RA and Dec arms, and the two polarizations, should agree within the
  errors (``CONS_ARM`` and ``CONS_POL`` are their differences in standard
  deviations; a warning is printed above ``--nsigma``). A large difference
  between arms can mean a pointing offset or a nearby source.
- ``CONS_SHP`` compares the default fit (baseline and beam fitted
  together) with a fit of the baseline outside a window around the source.
  A large negative value (a warning is printed) means that the source is
  wider than the beam or offset from the nominal position, and its flux is
  probably underestimated. A large positive value means that the beam is
  narrower than assumed (see ``--beam-fwhm``) or that there is structure
  near the source.
- The forced fit measures the peak of a point source at the nominal
  position. A pointing offset, or a source wider than the beam, make the
  flux lower. If the profiles are visibly wider than the fitted beam, or
  their peak is off-centre, the flux is a lower limit to the true one.
- Other sources along the scans, outside the region around the target, can
  bend the baseline and bias the flux (by up to half of their own flux).
  See :doc:`../technical_details` for the method and its validation.
- ``--source-coords`` measures the flux at a different position, e.g. to
  check that the result is zero away from the source.
