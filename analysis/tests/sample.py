"""Copy data files and keep a number of randomly sampled lines.

Used to create smaller test data from real data.

"""
import sys
import os
import itertools
import random

thisFileName = 'sample.py'

if not len(sys.argv) == 2:
    raise ValueError(
        f"{thisFileName} takes one argument, the top directory containing "
        "files to truncate."
    )

# Top directory of files to truncate.
topDirectory = sys.argv[1]
# Top directory to save truncated files to.
outDirectory = 'sampled'
# Number of rows to include at the start of the file.
includeRows = 2
# Number of lines of data to sample from each file.
numLines = 2000
# Directories to skip.
skipDirs = ['analysis', '.git']
# Files to skip.
skipFiles = ['.gitignore']

if os.path.isdir(outDirectory):
    raise OSError(
        f"{thisFileName} writes to the directory '{outDirectory}'"
        f" and won't overwrite it. Check if it exists already."
    )

for root, dirs, files in os.walk(topDirectory):
    # Skip files and directories given above.
    dirs[:] = [d for d in dirs if d not in skipDirs]
    files[:] = [f for f in files if f not in skipFiles]
    # Make output directory.
    rootRel = os.path.relpath(root, topDirectory)
    os.mkdir(os.path.join(outDirectory, rootRel))

    for file in files:
        outFile = os.path.join(outDirectory, rootRel, file)
        with open(os.path.join(root, file), 'r') as f:
            lines = f.readlines()
            header = lines[:includeRows]
            # Sample lines.
            sampled = random.sample(lines[includeRows:], numLines)
            with open(outFile, 'w') as outF:
                outF.write(''.join(header + sampled))
