"""Some nice helper functions."""

import os
import math
import numpy as np

def doublePrint(*objs, filePath=None, set=False, **kwargs):
    """Print a message and append it to a file at the same time!
    
     `objs`: Object(s) to print and write.
    `filePath`: File to write to. Doesn't write to a file by default.
    
    `end` and keyword arguments are the same as in the built-in
    `print`.
    
    """
    if set:
        doublePrint.filePath = filePath
    else:
        print(*objs, **kwargs)
        if filePath:
            message = ' '.join(objs) + '\n'
            with open(filePath, 'a') as file:
                file.write(message)
    # except Exception as e:
    #     print(f"Couldn't write to {filePath}!")
    #     print(e.message, e.args)
    #     print('Moving on anyway...')


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
    """Make getting the midpoints of an array of bin endpoints quicker."""
    return (bins[:-1] + bins[1:]) / 2
