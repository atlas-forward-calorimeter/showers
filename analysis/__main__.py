"""This runs when the analysis package is run from the command line."""

import sys
import os

# Add the analysis package to the system path.
sys.path.insert(
    0, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
)

from analysis import MultiRun

args = sys.argv[1:]  # Command line arguments.
if not args:
    # Default behavior with no arguments.
    # Currently runs a test.
    mr = MultiRun('analysis\\tests\\data')
    mr.go()
elif len(args) == 1:
    # Command line argument should be the path of the directory to analyze.
    mr = MultiRun(args[0])
    mr.go()
