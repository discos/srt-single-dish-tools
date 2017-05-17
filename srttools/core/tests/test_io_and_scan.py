# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function)

from ..read_config import read_config
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import pytest

from ..scan import Scan
from ..io import print_obs_info_fitszilla
from ..io import locations
import os
import numpy as np
import glob


class Test1_Scan(object):
    @classmethod
    def setup_class(klass):
        import os
        global DEBUG_MODE
        DEBUG_MODE = True

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, 'data')

        klass.fname = \
            os.path.abspath(
                os.path.join(klass.datadir, 'gauss_dec',
                             'Dec0.fits'))
        h5file = klass.fname.replace('.fits', '.hdf5')
        if os.path.exists(h5file):
            os.unlink(h5file)

        klass.config_file = \
            os.path.abspath(os.path.join(klass.datadir, 'test_config.ini'))
        print(klass.config_file)

        read_config(klass.config_file)

    def test_print_info(self, capsys):
        print_obs_info_fitszilla(self.fname)
        out, err = capsys.readouterr()
        assert 'bandwidth' in out.lower()

    def test_scan(self):
        '''Test that data are read.'''

        scan = Scan(self.fname)

        scan.write('scan.hdf5', overwrite=True)
        scan2 = Scan('scan.hdf5')
        assert scan.meta == scan2.meta

    def test_coordinate_conversion_works(self):
        scan = Scan(self.fname)
        obstimes = Time(scan['time'] * u.day, format='mjd', scale='utc')
        ref_coords = SkyCoord(ra=scan['ra'][:,0],
                              dec=scan['dec'][:,0],
                              obstime=obstimes,
                              location=locations[scan.meta['site']]
                              )
        altaz = ref_coords.altaz
        assert np.allclose(altaz.az.rad, np.array(scan['az'][:,0]))
        assert np.allclose(altaz.alt.rad, np.array(scan['el'][:,0]))



    @classmethod
    def teardown_class(klass):
        """Cleanup."""
        os.unlink('scan.hdf5')


class Test2_Scan(object):
    @classmethod
    def setup_class(klass):
        import os
        global DEBUG_MODE
        DEBUG_MODE = True

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, 'data')

        klass.fname = \
            os.path.abspath(
                os.path.join(klass.datadir, 'spectrum',
                             'roach_template.fits'))

        h5file = klass.fname.replace('.fits', '.hdf5')
        if os.path.exists(h5file):
            os.unlink(h5file)
        klass.config_file = \
            os.path.abspath(os.path.join(klass.datadir, 'spectrum.ini'))
        print(klass.config_file)

        read_config(klass.config_file)

    def test_scan(self):
        '''Test that data are read.'''

        scan = Scan(self.fname)

        scan.write('scan.hdf5', overwrite=True)
        scan.baseline_subtract('rough')

    @pytest.mark.parametrize('fname', ['srt_data.fits', 'med_data.fits'])
    def test_coordinate_conversion_works(self, fname):
        scan = Scan(os.path.join(self.datadir, 'spectrum', fname))
        obstimes = Time(scan['time'] * u.day, format='mjd', scale='utc')
        ref_coords = SkyCoord(ra=scan['ra'][:,0],
                              dec=scan['dec'][:,0],
                              obstime=obstimes,
                              location=locations[scan.meta['site']]
                              )
        altaz = ref_coords.altaz
        assert np.allclose(altaz.az.rad, np.array(scan['az'][:,0]))
        assert np.allclose(altaz.alt.rad, np.array(scan['el'][:,0]))

    @classmethod
    def teardown_class(klass):
        """Cleanup."""
        os.unlink('scan.hdf5')
        for f in glob.glob(os.path.join(klass.datadir, 'spectrum', '*.pdf')):
            os.unlink(f)
        for f in glob.glob(os.path.join(klass.datadir, 'spectrum', '*.hdf5')):
            os.unlink(f)

