"""Analyzing data in "pieces. Events, Runs, etc."""

import abc
import os

import pandas

from . import calc

_e_lims_200 = (0, 350)  # for 200 GeV
_e_lims_350 = (0, 700)  # for 350 GeV
_tube_e_lims_200 = tuple(lim / 15 for lim in _e_lims_200)
_tube_e_lims_350 = tuple(lim / 15 for lim in _e_lims_350)

# Skip this many lines at the start when reading hits data files.
_skiprows = 1


class Piece(abc.ABC):
    """
    A "piece" of analysis.

    Attributes:
        # TODO: Document this.
        hists: Contains all of the Histograms.
    """

    def __init__(self, input_path, out_dir=None, params=None):
        """

        :param input_path: Path to input file or directory.
        :type input_path: str
        :param out_dir: Path to output directory.
        :type out_dir: str or None
        :param params: Specific parameters for this run or event.
        :type params: dict or None
        """
        self._check_input_path(input_path)

        self.input_path = input_path
        self.out_dir = out_dir
        self.params = params

        self.name = self._input_path2name(input_path)
        self.param_values = self._parse_params(params)

        self.hists = {}
        # self.numbers = None

    @abc.abstractmethod
    def _check_input_path(self, input_path):
        """
        Make sure the input path is good.

        :param input_path:
        :type input_path: str
        :return:
        :rtype: None
        """

    @abc.abstractmethod
    def _parse_params(self, params):
        """
        Turn dictionary parameters into usable values.
        :param params: See __init__.
        :type params:
        :return: Usable values.
        :rtype:
        """

    @abc.abstractmethod
    def _input_path2name(self, input_path):
        """
        Turn an input path into a nicer name.

        :param input_path:
        :type input_path: str
        :return: A nicer name.
        :rtype: str
        """


class Event(Piece):
    """
    A single event.

        Attributes:
            name: A nice name to go by.

            hits: Data obtained from the hits file path.

            energy_vs_z: Energy vs. z histogram object, from calc.

            energy_vs_xy: See calc.

            numbers: See calc.

            run: The parent run that this event belongs to.
    """

    def __init__(self, hits_path, out_dir=None, params=None,
                 run=None):
        """

        :param hits_path: Path to particle hits data from Geant4.
        :type hits_path: str
        :param out_dir: Path to output directory.
        :type out_dir: str or None
        :param out_text: Path to output text. Defaults to
            out_dir/analysis.txt.
        :type out_text: str or None
        :param params: Specific parameters for this event.
        :type params: dict or None
        :param run: The parent run that this event belongs to.
        :type run: Run or None
        """
        super().__init__(hits_path, out_dir, params)

        self.run = run

        self.hits = None

        e_lims, tube_e_lims = self.param_values
        self.hists['energy_vs_z'] = calc.EnergyVsZ(
            e_lims, tube_e_lims, save_name=self.name, save_dir=self.out_dir
        )
        self.hists['energy_vs_xy'] = calc.EnergyVsXY(
            save_name='_'.join((self.name, 'xy')), save_dir=self.out_dir
        )
        self.numbers = calc.Numbers('numbers')

        self.go()

    def _check_input_path(self, input_path):
        assert os.path.isfile(input_path), \
            "Event data path isn't a file."

    def _parse_params(self, params):
        # TODO: Take out this test.
        if params is None:
            params = {'incident_energy': '200GeV'}

        if params['incident_energy'] == '200GeV':
            e_lims = _e_lims_200
            tube_e_lims = _tube_e_lims_200
        else:
            e_lims = _e_lims_350
            tube_e_lims = _tube_e_lims_350

        return e_lims, tube_e_lims

    def _input_path2name(self, input_path):
        return _file2name(input_path)

    def go(self):
        """Read in data, make plots and do calculations."""
        self.hits = pandas.read_csv(self.input_path, skiprows=_skiprows)

        for histogram in self.hists.values():
            histogram.add_data(self.hits)
            histogram.plot_single()

        tags = {'event': self.name}
        if self.run is not None:
            tags['run'] = self.run.name

        self.numbers.add_data(self.hits, tags)

    def __repr__(self):
        return f'<Event: {self.name} @ {self.input_path}>'


class Run(Piece):
    """A collection of events with a given setup."""

    def __init__(self, events_path, out_dir=None, params=None):
        """

        :param events_path: Path to directory containing data from
            multiple events.
        :type events_path: str
        :param out_dir: Save things in here.
        :type out_dir: str or None
        :param params: Specific parameters for this run.
        :type params: dict or None
        """
        super().__init__(events_path, out_dir, params)

        # TODO: Re-organize directory creation, and Run/MultiRun parent.
        if self.out_dir is None:
            self._out_subdir = None
        else:
            self._out_subdir = os.path.join(self.out_dir, self.name)
            os.makedirs(self._out_subdir, exist_ok=True)

        e_lims, tube_e_lims = self.param_values
        self.hists['energy_vs_z'] = calc.EnergyVsZ(
            e_lims, tube_e_lims, save_name=self.name, save_dir=self.out_dir
        )
        self.hists['energy_vs_xy'] = calc.EnergyVsXY(
            save_name='_'.join((self.name, 'xy')), save_dir=self.out_dir
        )
        self.numbers = calc.Numbers(save_name=self.name, save_dir=self.out_dir)

        self.go()

    def go(self):
        """Create new analyzed Events and analyze this Run."""
        self._new_events()

        for hist in self.hists.values():
            hist.plot_single(0)

        self.numbers.append_mean_and_dev(
            mean_tags={'run': self.name, 'event': 'mean'},
            dev_tags={'run': self.name, 'event': 'std_dev'}
        )
        self.numbers.save()

    def _new_events(self):
        """Create new analyzed Events and update this Run."""
        for hits_path in self._event_paths():
            event = Event(
                hits_path,
                out_dir=self._out_subdir,
                params=self.params,
                run=self
            )

            for hist_name, hist in event.hists.items():
                self.hists[hist_name].add_results(hist.get())

            self.numbers.add_results(event.numbers.get())

    def _event_paths(self):
        """Get the paths to the event data files."""
        root, dirs, files = next(os.walk(self.input_path))
        for filename in files:
            yield os.path.join(root, filename)

    def _check_input_path(self, input_path):
        assert os.path.isdir(input_path), \
            "Run directory path isn't actually a directory."

    def _parse_params(self, params):
        # TODO: Take out this test.
        if params is None:
            params = {'incident_energy': '200GeV'}

        if params['incident_energy'] == '200GeV':
            e_lims = _e_lims_200
            tube_e_lims = _tube_e_lims_200
        else:
            e_lims = _e_lims_350
            tube_e_lims = _tube_e_lims_350

        return e_lims, tube_e_lims

    def _input_path2name(self, input_path):
        return _dir2name(input_path)


def _file2name(file):
    """
    Turn a file path into a nicer name.

    :param file: File path.
    :type file: str
    :return: Nice name.
    :rtype: str
    """
    tail, head = os.path.split(file)
    assert head != '', "Is this a directory instead of a file?"

    return head.split('.')[0]


def _dir2name(directory):
    """
    Turn a directory path into a nicer name.

    :param directory: Directory path.
    :type directory: str
    :return: Nice name.
    :rtype: str
    """
    tail, head = os.path.split(directory)
    if head == '':
        tail, head = os.path.split(tail)

    return head
