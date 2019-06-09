"""The `Event` class."""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from .fcalpiece import FCalPiece
from .helpers import doublePrint, outTextPath

class Event(FCalPiece):
    """A collection of runs.
    
    One level higher than a run. The top level.
    
    """
    def __init__(self, filePath, outDirectory, params, name=None, parent=None):
        """ filePath: Path to this event's data file.
        outDirectory: Write all analysis output in here.

        For the docstrings of `params`, `maxEvents`, `name`, and 
        `parent`, see the constructor of `MultiRun`.

        """
        super().__init__(filePath, outDirectory, isAtom=True, parent=parent)

        self.filePath = filePath  # synonym for `inputPath`
        self.params = params

        self.fullEdep = 0
        self.middleEdep = 0
    
    def analyze(self):
        """Calculations and plots per event."""
        df = pd.read_csv(self.filePath, skiprows=1)

        # Energy deposit over all data.
        self.fullEdep = df.energy_deposit.sum()
        # Energy deposit over middle tube electrode.
        self.middleEdep = df.energy_deposit[
            df.z.abs() < self.params['tube middle z']
        ].sum()

        # Calculate energy-z histograms.
        self.fullSums, _ = np.histogram(df.z, 
                                   bins=self.parent.fullBins,
                                   weights=df.energy_deposit)
        # Calculate 2D histograms.
        dfMiddle = df[df.z.abs() < self.params['tube middle z']]
        self.xySums, _, _ = np.histogram2d(
            dfMiddle.x, 
            dfMiddle.y, 
            bins=(self.parent.xyBins, self.parent.xyBins), 
            weights=dfMiddle.energy_deposit
        )

        # Plot histograms.
        plt.plot(self.parent.fullBinMids, self.fullSums, lw=0.5)

        # Print stuff.
        doublePrint(self.filePath)
        doublePrint(f'fullEdep: {self.fullEdep}')
        doublePrint(f'middleEdep: {self.middleEdep}')
        