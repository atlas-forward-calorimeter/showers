"""Some nice helper functions."""

import os
import math
import numpy as np

def printAndWrite(*objs, **kwargs):
    """Print a message and append it to a file at the same time!
    
    Print once to the terminal, and then, if the `file` keyword argument
    is given, print to that file. If `file` is a string, it will be
    treated as a file path and opened in append mode. Otherwise, `file`
    is left as is when given to the built-in `print` function.
    
    """
    file = kwargs.pop('file', None)
    print(*objs, **kwargs)

    if file:
        if isinstance(file, str):
            with open(file, 'a') as f:
                print(*objs, file=f, **kwargs)
        else:
            print(*objs, file=file, **kwargs)


def outTextPath(outDirectory):
    if outDirectory:
        return os.path.join(outDirectory, 'analysis.txt')
    else:
        return None


def makeBins(start, end, binDensity=None):
    """Make histogram bins with a given bin density.

    Histogram bin density is fixed, so the bins begin exactly at `start` and 
    end at a rounded up value near `end`.
    """
    assert(end > start)

    if binDensity is None:
        binDensity = 10  # default

    # Round the end limit so the bin size is correct.
    numBins = math.ceil(binDensity * (end - start))
    newEnd = start + numBins / binDensity

    return np.linspace(start, newEnd, numBins + 1)


def binMidpoints(bins):
    """Get midpoints of bins from an array of bin endpoints."""
    return (bins[:-1] + bins[1:]) / 2
