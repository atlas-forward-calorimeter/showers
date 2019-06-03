"""Nice utility functions for use in the rest of the package.

Also provides access to some external packages (os, math, numpy, pandas).

Written by Anson Kost with the help of Professor John Rutherfoord, May 2019.

"""

import os
import math
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


def readFolder(folderPath, maxFiles=None, skiprows=1):
    """Read all of the events from a run."""
    if not maxFiles:  # Change 0, etc. to `None`.
        maxFiles = None
    
    folderPath = os.path.normpath(folderPath)  # Clean path format.
    
    _, subFolders, files = next(os.walk(folderPath))
    if files:
        # If files exist, don't read subfolders.
        for fileName in files[:maxFiles]:
            yield (
                fileName, 
                pd.read_csv(os.path.join(folderPath, fileName), 
                            skiprows=skiprows)
            )
    else:
        # Otherwise, keep going down.
        for folderName in subFolders:
            yield folderName, readFolder(os.path.join(folderPath, folderName))

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


def makeHist(df, bins=None, start=None, end=None, binDensity=None):
    """Make and save an energy-z histogram.
    
    If no bins are given, then bins are created using `makeBins`.
    In this case, if `start` and `end` are not given, they default to the min 
    and max of the dataframe's z values.
    """
    # Default bins.
    if bins is None:
        # Default limits cover all data.
        if start is None:
            start = df.z.min()
        if end is None:
            end = df.z.max()
        bins = makeBins(df.z.min(), df.z.max(), binDensity)

    sums, bins = np.histogram(df.z, bins=bins, weights=df.energy_deposit)
    '''
    if showPlot or savePath:
        binMidpoints = (bins[:-1] + bins[1:]) / 2
        plt.figure()
        plt.plot(binMidpoints, sums)
        if savePath:
            plt.savefig(savePath)
        if not showPlot:
            plt.close()
    '''
    return sums, bins
