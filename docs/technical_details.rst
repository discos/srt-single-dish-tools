Technical details
=================

This page collects implementation and maintenance notes that are useful to
developers, and the reasoning behind some non-obvious choices.

Supported Python and dependency versions
----------------------------------------

The package supports Python 3.9 and newer, because some production machines
at the telescope still run Python 3.9.

The lower bounds of the dependencies in ``pyproject.toml`` are not guesses:
they are the oldest versions that install from binary wheels on Python 3.9
and pass the test suite. They are exercised by the ``py39-test-oldestdeps``
tox environment, which uses ``uv`` with ``--resolution lowest-direct`` to
install exactly the lower bound of every direct dependency.

Some bounds are dictated by availability of binary wheels rather than by
features used in the code:

* ``scipy>=1.7.3``, ``h5py>=3.7``, ``pyyaml>=6.0``: older versions have no
  wheels for some current platforms (e.g. Apple Silicon) and would need to be
  compiled.
* ``matplotlib>=3.6``: ``rfistat`` uses ``width_ratios``/``height_ratios``
  as arguments of ``plt.subplots``.

Known issue: on Apple Silicon, the combination ``numpy==1.21.0`` +
``scipy==1.7.3`` can crash (segmentation fault inside the SVD used by
``curve_fit``) in the opacity tests. This is a problem of the bundled
OpenBLAS libraries of those old wheels, not of this package; the
oldest-dependency job in CI runs on Linux.

Running the tests
-----------------

The recommended way is through tox, with the ``tox-uv`` plugin:

.. code-block:: console

    $ pip install tox tox-uv
    $ tox -e py313-test-alldeps     # latest versions, all optional dependencies
    $ tox -e py39-test-oldestdeps   # oldest supported versions
    $ tox -e py313-test-devdeps     # nightly builds of numpy, scipy, matplotlib, astropy
    $ tox -e codestyle              # pre-commit checks (ruff, codespell, ...)

Tox installs the package in a separate environment and runs the tests from
``.tmp/<envname>``, so that the source tree is not imported by mistake. On
macOS, if tox picks up a "universal2" Python from python.org, the oldest
versions of numpy will not install; ``UV_PYTHON_PREFERENCE=only-managed``
forces tox-uv to use a uv-managed Python instead.

Warning policy in tests
-----------------------

``DeprecationWarning`` is turned into an error only when it is attributed to
``srttools`` itself (filter ``error::DeprecationWarning:srttools.*`` in
``pyproject.toml``). In practice, this catches calls from our code to
deprecated APIs of numpy, scipy, astropy etc., while ignoring deprecation
warnings raised *between* third-party packages (e.g. Pillow or pyparsing
warnings triggered inside old matplotlib versions), which we cannot fix.

Dates and time zones
--------------------

Dates given on the command line (e.g. ``SDTinspect --only-after
20250101-000000``) are always interpreted as UTC, independent of the time
zone of the machine. Tests force a non-UTC time zone to check this, because
CI machines run in UTC and would not catch the problem otherwise.
