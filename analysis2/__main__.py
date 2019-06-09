import sys
import os

# Fix path so relative imports from this package are okay.
# sys.path.insert(
#     0, 
#     os.path.abspath(os.path.dirname(__file__))
# )
## Above doesn't work; take it out if there's no progress.

from core.multirun import MultiRun

args = sys.argv[1:]  # Command line arguments.
if not args:
    # Default behavior with no arguments.
    # Currently runs a test.
    mr = MultiRun('analysis2\\tests\\test_data2')
    mr.go()
else:
    # First CL argument should be the path of the directory to analyze.
    pass
