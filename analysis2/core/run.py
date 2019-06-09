"""The `Run` class."""

import os
import numpy as np
from matplotlib import pyplot as plt

from .fcalpiece import FCalPiece
from .event import Event
from .helpers import doublePrint, outTextPath, makeBins, binMidpoints

class Run(FCalPiece):
    """A collection of runs.
    
    One level higher than a run. The top level.
    
    """
    def __init__(self, 
                 dataDirectory, 
                 outDirectory,
                 params,
                 maxEvents=None,
                 name=None, 
                 parent=None):
        """dataDirectory: Directory of the run's data. Should contain 
                      event data files.
        outDirectory: Write all analysis output in here.
        
        For the docstrings of `params`, `maxEvents`, `name`, and 
        `parent`, see the constructor of `MultiRun`.

        """
        super().__init__(dataDirectory, outDirectory, parent=parent)

        self.dataDirectory = dataDirectory  # synonym for `inputPath`
        self.params = params
        self.maxEvents = maxEvents

        ## Initialize analysis stuff. ##

        # Histogram bins.
        self.fullBins = makeBins(*self.params['full z limits'],
                                 self.params['bin density'])
        self.fullBinMids = binMidpoints(self.fullBins)
        self.xyBins = makeBins(-self.params['plate xy'],
                               self.params['plate xy'],
                               self.params['bin density'])
        self.xyBinMids = binMidpoints(self.xyBins)

        # Histogram sums.
        self.fullSums = np.zeros(len(self.fullBinMids))
        self.xySums = np.zeros((len(self.xyBinMids), 
                                   len(self.xyBinMids)))

        # Energy sums.
        self.fullEdeps, self.middleEdeps = [], []

        self.intialize()
    
    def intialize(self):
        """Like a constructor, but for analysis and output."""        
        # Run Header
        runHeader = f'\nBegin of run {self.dataDirectory}.\n'
        self.doublePrint(runHeader)

        # Multiple event histogram figure and formatting.
        plt.figure()
        plt.title('Histogram - Energy Deposit vs. z'
                  f' - Events - Run {self.name}')
        plt.xlabel('z')
        plt.ylabel('Energy Deposit Per Bin')
    
    def analyzeSmaller(self, event):
        """Average over events."""
        self.fullEdeps.append(event.fullEdep)
        # self.fullSums += event.fullSums
        # self.xySums += event.xySums

    def analyze(self):
        """Finish up and output run calculations and plots."""
   
    def smallerPieces(self):
        """Get every event file in `dataDirectory`."""
        root, dirs, files = next(os.walk(self.dataDirectory))
        for file in files[:self.maxEvents]:  # Max number of events.
            yield Event(filePath=os.path.join(root, file),
                        outDirectory=self.outDirectory,
                        params=self.params,
                        parent=self)
