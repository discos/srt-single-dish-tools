import os
import glob
import pytest
from srttools.imager import main_preprocess
from srttools.rfistat import main_rfistat
from srttools.utils import HAS_MPL


class TestRFI:
    @classmethod
    def setup_class(klass):
        klass.datadir = os.path.join(os.path.dirname(__file__), "data")

    def test_rfi_config_mpl(self):
        config_file = os.path.abspath(os.path.join(self.datadir, "nodding_xarcos.ini"))

        main_preprocess(["-c", config_file])
        for f in glob.glob("rfi*jpg") + glob.glob("*rfi.hdf5"):
            os.unlink(f)
        if HAS_MPL:
            main_rfistat(["-c", config_file])
            assert len(glob.glob("rfi*jpg")) > 0
        else:
            with pytest.warns(UserWarning, match="Matplotlib is required"):
                main_rfistat(["-c", config_file])

        assert len(glob.glob("*rfi.hdf5")) > 0
        for f in glob.glob("rfi*jpg") + glob.glob("*rfi.hdf5"):
            os.unlink(f)

    def test_rfi_files_mpl(self):
        files = glob.glob(os.path.join(os.path.join(self.datadir, "nodding_xarcos"), "*.hdf5"))
        if HAS_MPL:
            main_rfistat(files)
            assert len(glob.glob("rfi*jpg")) > 0
        else:
            with pytest.warns(UserWarning, match="Matplotlib is required"):
                main_rfistat(files)

        assert len(glob.glob("*rfi.hdf5")) > 0
        for f in glob.glob("rfi*jpg") + glob.glob("*rfi.hdf5"):
            os.unlink(f)
