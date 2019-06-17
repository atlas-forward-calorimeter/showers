"""The `Event` class."""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from .piece import FCalPiece
from .helpers import printAndWrite

class Event(FCalPiece):
    """A collection of runs.
    
    One level higher than a run. The top level.
    
    """
    def __init__(self, filePath, outDirectory=None, parent=None):
        """filePath: Path to this event's data file."""
        super().__init__(filePath, outDirectory=outDirectory, parent=parent)

        self.filePath = filePath  # same as `self.inputPath`

        self.fullEdep = 0
        self.middleEdep = 0
        self.zFullSums = None
        self.xyMiddleSums = None
    
    def analyze(self):
        """Calculations and plots per event."""
        df = pd.read_csv(self.filePath, skiprows=1)
        dfMiddle = df[df.z.abs() < self.tubeMiddleZ]

        # Energy deposit over all data.
        self.fullEdep = df.energy_deposit.sum()
        # Energy deposit over middle tube electrode.
        self.middleEdep = dfMiddle.energy_deposit.sum()

        # Calculate energy-z histograms.
        self.zFullSums, _ = np.histogram(
            df.z, bins=self.fullBins, weights=df.energy_deposit)
        # Calculate 2D histograms.
        self.xyMiddleSums, _, _ = np.histogram2d(
            dfMiddle.x, 
            dfMiddle.y, 
            bins=2 * (self.xyBins,), 
            weights=dfMiddle.energy_deposit)
        
        # Update the run.
        if self.parent:
            self.parent.analyzeEvent(self)

        # Plot event histogram.
        plt.plot(self.fullBinMids, self.zFullSums, lw=0.5)

        # Print stuff.

        output = (
            f'Event {self.name}.\n'
            f'fullEdep: {self.fullEdep}.\n'
            f'middleEdep: {self.middleEdep}.'
        )
        printAndWrite(output, file=self.outTextPath)
