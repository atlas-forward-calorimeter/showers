"""The `Run` class."""

import os
import numpy as np
from matplotlib import pyplot as plt

from .piece import FCalPiece
from .event import Event
from .helpers import printAndWrite, makeBins, binMidpoints

class Run(FCalPiece):
    """A run, which is a collection of events with a given setup.

    The important graphs and calculations will be of averages per run.
    
    """
    def __init__(self, 
                 dataDirectory, 
                 outDirectory=None,
                 maxEvents=None,
                 parent=None):
        """dataDirectory: Directory of the run's data. Should contain 
                      event data files.
        
        For the docstrings of `maxEvents` and `parent`, see the 
        constructor of `MultiRun`.

        """
        super().__init__(
            dataDirectory, 
            outDirectory=outDirectory, 
            maxEvents=maxEvents, 
            parent=parent)

        self.dataDirectory = dataDirectory  # same as `self.inputPath`

        ## Initialize analysis stuff. ##

        # Histogram sums.
        self.zFullSums = np.zeros(len(self.fullBinMids))
        self.xyMiddleSums = np.zeros(
            2 * (len(self.xyBinMids),)
        )

        # Energy sums.
        self.fullEdeps, self.middleEdeps = [], []

        # Results.
        self.meanEdep = 0
        self.eDepSigma = 0
        self.meanMiddleEdep = 0
        self.middleEdepSigma = 0

        self.start()
    
    def start(self):
        """Like a constructor, but for analysis and output."""        
        # Run Header
        runHeader = f'\nBegin of run {self.name}.\n'
        printAndWrite(runHeader, file=self.outTextPath)

        # Multiple event histogram figure and formatting.
        plt.figure()
        plt.title('Histogram - Energy Deposit vs. z'
                  f' - Events - Run {self.name}')
        plt.xlabel('z')
        plt.ylabel('Energy Deposit Per Bin')
    
    def analyzeEvent(self, event):
        """Update the run per event."""
        self.fullEdeps.append(event.fullEdep)
        self.middleEdeps.append(event.middleEdep)
        self.zFullSums += event.zFullSums
        self.xyMiddleSums += event.xyMiddleSums

    def analyze(self):
        """Calculate averages over the run and plot them."""
        numEvents = len(self.fullEdeps)
        fullEdeps = np.array(self.fullEdeps)
        middleEdeps = np.array(self.middleEdeps)

        self.meanEdep = np.mean(fullEdeps)
        self.eDepSigma = np.std(fullEdeps)
        self.meanMiddleEdep = np.mean(middleEdeps)
        self.middleEdepSigma = np.std(middleEdeps)

        # Run calculations output.
        output = (
            '-\n'
            f'Run {self.name}.\n'
            'Total Energy Deposit Average:\n'
            f'{self.meanEdep} ({self.eDepSigma}) MeV\n'
            'Middle Tube Energy Deposit Average:\n'
            f'{self.meanMiddleEdep} ({self.middleEdepSigma}) MeV\n'
            '-\n'
        )

        # Run Footer
        footer = f'\nEnd of run {self.name}.\n'
        printAndWrite(output + footer, file=self.outTextPath)

        # Finish run histograms.

        # Event energy-z.
        eventHistFilename = f'{self.name}-EventHist.{self.plotFileFormat}'
        if self.outDirectory:
            plt.savefig(
                os.path.join(self.outDirectory, eventHistFilename), 
                format=self.plotFileFormat)
        # Energy-z.
        plt.figure()
        plt.title(f'Histogram - Energy Deposit vs. z - Sum - Run {self.name}')
        plt.xlabel('z')
        plt.ylabel('Energy Deposit Per Bin')
        plt.plot(self.fullBinMids, self.zFullSums)
        ezHistFilename = f'{self.name}-Hist.{self.plotFileFormat}'
        if self.outDirectory:
            plt.savefig(
                os.path.join(self.outDirectory, ezHistFilename), 
                format=self.plotFileFormat)
        # Energy-xy.
        plt.figure()
        plt.title(f'Histogram - Energy Deposit vs. x, y - Sum - Run {self.name}')
        plt.xlabel('x')
        plt.ylabel('y')
        gridX, gridY = np.meshgrid(self.xyBins, self.xyBins)
        plt.pcolormesh(gridX, gridY, self.xyMiddleSums)
        # if self.outDirectory:
        #    xyHistFilename = f'{self.name}-xyHist.{self.plotFileFormat}'
        #    plt.savefig(
        #        os.path.join(self.outDirectory, xyHistFilename), 
        #        format=self.plotFileFormat)
   
    def smallerPieces(self):
        """Get every event file in `dataDirectory`."""
        root, dirs, files = next(os.walk(self.dataDirectory))
        for file in files[:self.maxEvents]:  # Max number of events.
            yield Event(filePath=os.path.join(root, file), parent=self)
