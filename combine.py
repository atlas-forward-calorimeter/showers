"""Combine existing results output into one nice table, from the command
line.

Usage:
    $ python combine.py [out_name] *[result_paths]

The means and uncertainties of the results from each run are compiled
into a single table. Columns for the `info` of the runs are also
created. Compiled results are saved to files named `out_name`. If no
`result_paths` are given, combine.py will try to combine all csv files
within the current directory and subdirectories.
"""

import os
import sys

import pandas

from analysis import pieces

# Skip this many lines at the start when reading analysis output files.
_skiprows = 1


def _csv_files(path):
    """
    Iterate over all csv files inside `path` and its subdirectories.
    """
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith('.csv'):
                yield file


if __name__ == '__main__':
    assert len(sys.argv) >= 2, \
        "Must provide an output file name ('none' for no output)."

    out_file = None if sys.argv[1].lower() == 'none' else sys.argv[1]

    if len(sys.argv) > 2:
        result_paths = sys.argv[2:]
    else:
        result_paths = _csv_files(os.getcwd())

    pieces.combine_results(
        (
            pandas.read_csv(file, skiprows=_skiprows, delim_whitespace=True)
            for file in result_paths
        ),
        out_file=out_file
    )
