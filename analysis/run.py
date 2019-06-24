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

        # Multiple event plot.
        self.fig = None
        self.ax = None

        # Histogram sums.
        self.zFullSumss = []
        self.xyMiddleSumss = []
        self.zMeanFullSums = None
        self.xyMeanMiddleSums = None
        self.zFullSumSigmas = None

        # Energy sums.
        self.fullEdeps, self.middleEdeps = [], []

        # Results.
        self.meanEdep = 0
        self.eDepSigma = 0
        self.meanMiddleEdep = 0
        self.middleEdepSigma = 0

        # Histogram limits.

        if "350GeV" in self.name:
            self.histYlim = 70
        else:
            self.histYlim = 35

        self.middleHistYlim = self.histYlim / 20

        self.start()
    
    def start(self):
        """Like a constructor, but for analysis and output."""
        if self.outDirectory:
            # Create output directory
            # (and overwrite any existing analysis).
            os.makedirs(self.outDirectory, exist_ok=True)

        # Run Header
        runHeader = f'\nBegin of run {self.name}.\n'
        printAndWrite(runHeader, file=self.parent.outTextPath)

        # Multiple event histogram figure and formatting.
        self.fig, self.ax = plt.subplots()
        self.ax.set_title('Histogram - Energy Deposit vs. z'
                  f' - Events - Run {self.name}')
        self.ax.set_xlabel('z')
        self.ax.set_ylabel('Energy Deposit Per Bin')
        self.ax.set_ylim(0, self.histYlim)
    
    def analyzeEvent(self, event):
        """Update the run per event."""
        self.fullEdeps.append(event.fullEdep)
        self.middleEdeps.append(event.middleEdep)
        self.zFullSumss.append(event.zFullSums)
        self.xyMiddleSumss.append(event.xyMiddleSums)

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
        printAndWrite(output + footer, file=self.parent.outTextPath)

        # Finish run histograms.

        zFullSumss = np.array(self.zFullSumss)
        xyMiddleSumss = np.array(self.xyMiddleSumss)
        # Bin sum averages.
        self.zMeanFullSums = np.mean(zFullSumss, axis=0)
        self.xyMeanMiddleSums = np.mean(xyMiddleSumss, axis=0)
        # Bin standard deviations.
        self.zFullSumsSigmas = np.std(zFullSumss, axis=0)
        # self.xyMiddleSumsSigmas = np.std(xyMiddleSumss, axis=0)

        # Multiple event energy-z.
        eventHistFilename = f'{self.name}-EventHist.{self.plotFileFormat}'
        if self.parent.outDirectory:
            self.fig.savefig(
                os.path.join(self.parent.outDirectory, eventHistFilename),
                format=self.plotFileFormat)

        # Close multiple event energy-z.
        plt.close(self.fig)

        # Energy-z.
        plt.figure()
        plt.title(f'Histogram - Energy Deposit vs. z - Sum - Run {self.name}')
        plt.xlabel('z')
        plt.ylabel('Energy Deposit Per Bin')
        plt.ylim(0, self.histYlim)
        plt.plot(self.fullBinMids, self.zMeanFullSums, lw=0.6)
        plt.plot(self.fullBinMids, self.zMeanFullSums - self.zFullSumsSigmas, lw=0.5)
        plt.plot(self.fullBinMids, self.zMeanFullSums + self.zFullSumsSigmas, lw=0.5)
        ezHistFilename = f'{self.name}-Hist.{self.plotFileFormat}'
        if self.parent.outDirectory:
            plt.savefig(
                os.path.join(self.parent.outDirectory, ezHistFilename),
                format=self.plotFileFormat)
        # Energy-xy.
        plt.figure()
        plt.title(f'Histogram - Energy Deposit vs. x, y - Sum - Run {self.name}')
        plt.xlabel('x')
        plt.ylabel('y')
        gridX, gridY = np.meshgrid(self.xyBins, self.xyBins)
        plt.pcolormesh(gridX, gridY, self.xyMeanMiddleSums)
        if self.parent.outDirectory:
           xyHistFilename = f'{self.name}-xyHist.{self.plotFileFormat}'
           plt.savefig(
               os.path.join(self.parent.outDirectory, xyHistFilename),
               format=self.plotFileFormat)
   
    def smallerPieces(self):
        """Get every event file in `dataDirectory`."""
        root, dirs, files = next(os.walk(self.dataDirectory))
        for file in files[:self.maxEvents]:  # Max number of events.
            yield Event(filePath=os.path.join(root, file), parent=self)
