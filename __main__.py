"""Handle arguments from the command line.

This module is run when the package itself is executed in the command
line, e.g., with the command:
    $ python ./showers/
"""

import sys
from analysis import pieces

if len(sys.argv) == 1:
    run = pieces.Run('.')
    print(run.numbers)
elif len(sys.argv) == 2:
    run = pieces.Run(sys.argv[1])
    print(run.numbers)
else:
    out_dir = None if sys.argv[1] == 'none' else sys.argv[1]
    runs = pieces.go(out_dir, *sys.argv[2:])
    for run in runs:
        print(run.numbers)
