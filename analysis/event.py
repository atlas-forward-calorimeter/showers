"""The Event class."""

import os
import numpy as np
import pandas
from matplotlib import pyplot as plt

from .piece import FCalPiece
from .helpers import printAndWrite
from . import calc

_e_lims_200 = (0, 350)  # for 200 GeV
_e_lims_350 = (0, 700)  # for 350 GeV
_tube_e_lims_200 = tuple(lim / 15 for lim in _e_lims_200)
_tube_e_lims_350 = tuple(lim / 15 for lim in _e_lims_350)

# Skip this many lines at the start when reading hits data files.
_skiprows = 1


class Event(FCalPiece):
    """A single event."""

    def __init__(self, hits_path, out_dir=None, out_text=None, params=None):
        """

        :param hits_path: Path to particle hits data from Geant4.
        :type hits_path: str
        :param out_dir: Path to output directory.
        :type out_dir: str
        :param out_text: Path to output text. Defaults to
            out_dir/analysis.txt.
        :type out_text: str
        :param params: Specific parameters for this event.
        :type params: dict

        Attributes:
            name: A nice name to go by.

            hits: Data obtained from the hits file path.

            energy_vs_z: Energy vs. z histogram object, from calc.

            energy_vs_xy: See calc.

            split_z: See calc.

            e_dep: Total energy deposit for all of this event's data.

            middle_e_dep: Energy deposit in all four middle tube
            sections.
        """
        assert os.path.isfile(hits_path)

        self._hits_path = hits_path
        self._out_dir = out_dir

        # public attributes
        self.name = _file2name(self._hits_path)
        self.hits = None
        self.energy_vs_z = None
        self.energy_vs_xy = None
        self.split_z = calc.SplitZ()
        self.e_dep = 0
        self.middle_e_dep = 0

        # defaults
        if (out_text is None) and (self._out_dir is not None):
            self._out_text = os.path.join(self._out_dir, 'analysis.txt')
        else:
            self._out_text = out_text

        if params is None:
            self._params = {'incident_energy': '350GeV'}
        else:
            self._params = params
        ####

        e_lims, tube_e_lims = _parse_params(self._params)
        self.energy_vs_z = calc.EnergyVsZ(
            e_lims, tube_e_lims, save_name=self.name, save_dir=self._out_dir
        )

        self.energy_vs_xy = calc.EnergyVsXY(
            save_name='_'.join((self.name, 'xy')), save_dir=self._out_dir
        )

        self.go()

        # super().__init__(hits_path, outDirectory=out_dir, parent=params)
        # self.filePath = hits_path  # same as `self.inputPath`
        """
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

        self.middleHistYlim = self.histYlim / 15

        self.start()
        """

    def go(self):
        # Read data.
        self.hits = pandas.read_csv(self._hits_path, skiprows=_skiprows)

        # Plot data.
        self.energy_vs_z.add_data(self.hits)
        self.energy_vs_xy.add_data(self.hits)
        self.energy_vs_z.plot_single()
        self.energy_vs_xy.plot_single()

        # More calcs.
        self.split_z.add_data(self.hits)

        # test
        print(self.split_z.get())

    def start(self):
        """Like a constructor, but for analysis and output."""
        self.fig, self.ax = plt.subplots()
        self.ax.set_title('Energy Deposit vs. z'
                          f'-Run {self.parent.name}-Event {self.name}')
        self.ax.set_xlabel(f'z ({self.lengthUnits})')
        self.ax.set_ylabel(f'Energy Deposit Per Bin ({self.energyUnits})')
        self.ax.set_ylim(0, self.histYlim)

        self.ax2 = self.ax.twinx()
        self.ax2.set_ylim(0, self.middleHistYlim)

    def analyze(self):
        """Calculations and plots per event."""
        df = pd.read_csv(self.filePath, skiprows=1)
        dfMiddle = df[df.z.abs() < self.tubeMiddleZ]

        # Energy deposit over all data.
        self.e_dep = df.energy_deposit.sum()
        # Energy deposit over middle tube electrode.
        self.middle_e_dep = dfMiddle.energy_deposit.sum()

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
            f'fullEdep: {self.e_dep}.\n'
            f'middleEdep: {self.middle_e_dep}.'
        )
        printAndWrite(output, file=self.parent.parent.outTextPath)

    def __repr__(self):
        return f'<Event: {self.name} @ {self._hits_path}>'


def _parse_params(params):
    """
    Parse parameters and return useful results.

    :param params: The parameters
    :type params: dict
    :return: Useful results.
    :rtype: any
    """
    if params['incident_energy'] == '200GeV':
        e_lims = _e_lims_200
        tube_e_lims = _tube_e_lims_200
    else:
        e_lims = _e_lims_350
        tube_e_lims = _tube_e_lims_350

    return e_lims, tube_e_lims


def _file2name(file):
    """

    :param file: File path.
    :type file: str
    :return: Nice name.
    :rtype: str
    """
    tail, head = os.path.split(file)
    assert head != '', "Is this a directory instead of a file?"

    return head.split('.')[0]
