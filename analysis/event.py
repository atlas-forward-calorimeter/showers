"""The `Event` class."""

import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from .piece import FCalPiece
from .helpers import printAndWrite


class Event(FCalPiece):
    """A single event."""

    def __init__(self, hits_data, out_dir=None, out_text=None, prefix=None):
        """

        :param hits_data: Path to particle hits data from Geant4.
        :type hits_data: str
        :param out_dir: Path to output directory.
        :type out_dir: str
        :param out_text: Path to output text.
        :type out_text: str
        :param prefix: Identifies the run of this event.
        :type prefix: str
        """
        assert os.path.isfile(hits_data)

        self.hits_data = hits_data
        self.out_dir = out_dir
        self.out_text = out_text
        self.prefix = prefix

        super().__init__(hits_data, outDirectory=out_dir, parent=prefix)

        self.filePath = hits_data  # same as `self.inputPath`

        self.fullEdep = 0
        self.middleEdep = 0
        self.zFullSums = None
        self.xyMiddleSums = None

        # Single event plot.
        self.fig = None
        self.ax = None
        self.ax2 = None

        # Histogram limits.

        if "350GeV" in self.parent.name:
            self.histYlim = 70
        else:
            self.histYlim = 35

        self.middleHistYlim = self.histYlim / 20

        self.start()

    def start(self):
        """Like a constructor, but for analysis and output."""
        self.fig, self.ax = plt.subplots()
        self.ax.set_title('Energy Deposit vs. z'
                          f'-Run {self.parent.name}-Event {self.name}')
        self.ax.set_xlabel('z')
        self.ax.set_ylabel('Energy Deposit Per Bin')
        self.ax.set_ylim(0, self.histYlim)

        self.ax2 = self.ax.twinx()
        self.ax2.set_ylim(0, self.middleHistYlim)

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
        tubeSlice = np.abs(self.fullBinMids) < self.tubeZ
        platesSlice = np.logical_not(tubeSlice)
        tubesPlot = (self.fullBinMids[tubeSlice], self.zFullSums[tubeSlice])
        platesPlot = (self.fullBinMids[platesSlice], self.zFullSums[platesSlice])

        self.parent.ax.plot(self.fullBinMids, self.zFullSums, lw=0.5)

        self.ax.plot(*platesPlot, lw=0.5)
        self.ax2.plot(*tubesPlot, lw=0.5)

        # Save it.
        singleEventHistFilename = f'{self.name}-Hist.{self.plotFileFormat}'
        if self.parent.outDirectory:
            self.fig.savefig(
                os.path.join(self.parent.outDirectory, singleEventHistFilename),
                format=self.plotFileFormat)

        # Close figure.
        plt.close(self.fig)

        # Print stuff.
        output = (
            f'Event {self.name}.\n'
            f'fullEdep: {self.fullEdep}.\n'
            f'middleEdep: {self.middleEdep}.'
        )
        printAndWrite(output, file=self.parent.parent.outTextPath)

        @staticmethod
        def filepath2name(filepath):
            head, tail = os.path.split(filepath)
            if not tail:
                _, tail = os.path.split(head)

