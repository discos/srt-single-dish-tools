from astropy.io import fits
import os
import numpy as np
from srttools.converters.mbfits import MBFITS_creator


class MBFits_template(object):
    def __init__(self, mbfits_dir='.'):
        self.mbfits_dir = mbfits_dir
        self.GROUPING = fits.open(os.path.join(mbfits_dir, 'GROUPING.fits'))
        self.files = self.GROUPING[1].data['MEMBER_LOCATION']
        self.extnames = self.GROUPING[1].data['EXTNAME']
        self.febes = self.GROUPING[1].data['FEBE']
        self.subsnums = self.GROUPING[1].data['SUBSNUM']
        self.basebands = self.GROUPING[1].data['BASEBAND']

        scan_file = self.files[self.extnames == 'SCAN-MBFITS'][0]
        self.SCAN = fits.open(os.path.join(mbfits_dir, scan_file))

        FEBE_combos = self.SCAN['SCAN-MBFITS'].data['FEBE']
        print(FEBE_combos)
        self.FEBEPARs ={}
        for febe in FEBE_combos:
            self.FEBEPARs[febe] = \
                fits.open(os.path.join(mbfits_dir, febe + '-FEBEPAR.fits'))
        self.time = None

    def read_subscan(self, file):
        hdul = fits.open(os.path.join(self.mbfits_dir, file))
        if self.time is None:
            self.time = hdul[1].data['MJD']
        else:
            try:
                if not np.allclose(self.time, hdul[1].data['MJD']):
                    raise ValueError('MJD mismatch in files')
            except ValueError:
                return None

        return hdul[1].data['DATA']

    def modify_subscan(self):
        pass

    def list_scans(self, febe, baseband):
        good = (self.febes == febe) & (self.basebands == baseband)
        return self.files[good]


if __name__ == '__main__':  # pragma: no cover
    import sys
    file = MBFits_template(sys.argv[1])
    febe = 'FLASH460L-XFFTS'
    print(file.list_scans(febe, 1))
    print(file.list_scans(febe, 2))
    files = file.list_scans(febe, 1)

    print(file.read_subscan(files[0]))

    created_file = MBFITS_creator('try')

    created_file.add_subscan(sys.argv[2])

    file = MBFits_template('try')

    febe = 'CCB0RCP-ROACH2'
    print(file.list_scans(febe, 1))
    print(file.list_scans(febe, 2))
    files = file.list_scans(febe, 1)

    print(file.read_subscan(files[0]))
