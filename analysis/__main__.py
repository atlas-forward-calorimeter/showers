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
else:
    # Command line arguments are passed to `MultiRun` constructor.
    mr = MultiRun(*args)
    mr.go()
