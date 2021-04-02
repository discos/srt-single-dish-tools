from srttools.parse_acs_logs import main_parse_acs_logs


class Test1_Scan(object):
    @classmethod
    def setup_class(klass):
        import os

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, "data")

        klass.fname = os.path.abspath(os.path.join(klass.datadir, "acs.xml"))

    def test_read_acs_log_list_calon(self, capsys):
        main_parse_acs_logs([self.fname, "--list-calon"])
        out, err = capsys.readouterr()
        assert "20180731-035506-Rescicom-N1068cal_001_028.fits" in out
        # There are four files in the output
        assert len([k for k in out.split("\n") if "2018" in k]) == 4
