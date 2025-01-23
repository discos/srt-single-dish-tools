import os
import glob
from srttools.imager import main_preprocess
from srttools.rfistat import main_rfistat


class TestRFI:
    @classmethod
    def setup_class(klass):
        klass.datadir = os.path.join(os.path.dirname(__file__), "data")

    def test_rfi_config(self):
        config_file = os.path.abspath(os.path.join(self.datadir, "nodding_xarcos.ini"))

        main_preprocess(["-c", config_file])
        for f in glob.glob("rfi*jpg") + glob.glob("*rfi.hdf5"):
            os.unlink(f)
        main_rfistat(["-c", config_file])
        assert len(glob.glob("rfi*jpg")) > 0
        assert len(glob.glob("*rfi.hdf5")) > 0
        for f in glob.glob("rfi*jpg") + glob.glob("*rfi.hdf5"):
            os.unlink(f)

    def test_rfi_files(self):
        files = glob.glob(os.path.join(os.path.join(self.datadir, "nodding_xarcos"), "*.hdf5"))
        main_rfistat(files)
        assert len(glob.glob("rfi*jpg")) > 0
        assert len(glob.glob("*rfi.hdf5")) > 0
        for f in glob.glob("rfi*jpg") + glob.glob("*rfi.hdf5"):
            os.unlink(f)
