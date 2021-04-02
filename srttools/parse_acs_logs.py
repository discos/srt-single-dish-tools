import xml.etree.ElementTree as ET
from astropy.table import Table
from astropy.time import Time
import os
import numpy as np


def load_acs_log_file(fname, full=False):
    with open(fname) as fd:
        doc = fd.read()
        # fix missing tags
        doc = "<!DOCTYPE xml>\n<html>\n" + doc + "\n</html>"
    tree = ET.fromstring(doc)

    is_cal_on = False

    table = Table(
        names=["kind", "Time", "file", "CAL", "text"],
        dtype=["<U10", "O", "<U200", bool, "<U200"],
    )

    for line_el in tree.iter():
        line = dict(line_el.attrib)
        if line == {}:
            continue
        if "TimeStamp" in line:
            time = Time(line["TimeStamp"])
        else:
            continue
        text = line_el.text
        if text is None:
            continue
        if "calOff" in text:
            is_cal_on = False
        elif "calOn" in text:
            is_cal_on = True

        file = ""
        if ".fits" in line_el.text:
            file = os.path.basename(text.replace("FILE_OPENED", "").strip())
            text = ""
        else:
            if not full:
                continue

        table.add_row([line_el.tag, time, file, is_cal_on, text])

    return table


def main_parse_acs_logs(args=None):
    import argparse

    description = "Read ACS logs and return useful information"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "files",
        nargs="*",
        help="Single files to preprocess",
        default=None,
        type=str,
    )
    parser.add_argument(
        "--to-csv",
        action="store_true",
        default=False,
        help="Save a CSV file with the results",
    )
    parser.add_argument(
        "--list-calon",
        action="store_true",
        default=False,
        help="List files with calibration mark on",
    )

    args = parser.parse_args(args)
    for fname in args.files:
        # Get the full table only when relevant
        full = args.to_csv
        table = load_acs_log_file(fname, full=full)
        if args.list_calon:
            print("\n\nList of files with the calibration mark on:\n")
            good = table["CAL"] & (table["file"] != "")
            if np.any(good):
                for fname in table["file"][good]:
                    print(fname)
