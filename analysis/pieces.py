"""Analyzing data in "pieces." Events, Runs, etc."""

import abc
import collections
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
        name: A nice name to go by.
        hists: Contains all of the Histograms.

    TODO: Document this.
    TODO: Parse paths and make titles.
    """

    def __init__(self, input_path, parent=None, out_dir=None, info=None):
        """

        :param input_path: Path to input file or directory.
        :type input_path: str
        :param out_dir: Path to output directory.
        :type out_dir: str or None
        :param info: Specific info for this run or event,
            e.g. incident energy.
        :type info: dict or None
        """
        self._check_input_path(input_path)

        self.input_path = input_path
        self.parent = parent
        self.out_dir = out_dir

        self.info = info or {}

        self.name = self._input_path2name(input_path)

        if self.parent:
            # Inherit things from the parent.
            if self.parent.info:
                _weak_update(self.info, self.parent.info)
            if not self.out_dir:
                self.out_dir = self.parent.out_subdir

        if self.out_dir:
            # Create out_dir in the filesystem and make a subdirectory
            # name.
            os.makedirs(self.out_dir, exist_ok=True)
            self.out_subdir = os.path.join(self.out_dir, self.name)
        else:
            self.out_subdir = None

        self.hists = {}

        e_lims, tube_e_lims = self._get_e_lims(self.info)
        self.hists['energy_vs_z'] = calc.EnergyVsZ(self, e_lims, tube_e_lims)
        self.hists['energy_vs_xy'] = calc.EnergyVsXY(self)

        self.numbers = calc.Numbers(self)

    def __repr__(self):
        return f'<{self.__class__.__name__ }: {self.name} @ {self.input_path}>'

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
    def _input_path2name(self, input_path):
        """
        Turn an input path into a nicer name.

        :param input_path:
        :type input_path: str
        :return: A nicer name.
        :rtype: str
        """

    @staticmethod
    def _get_e_lims(info):
        """
        Get energy plot limits (full and tubes) from Piece info.

        :param info: The Piece info.
        :return: Energy plot limits (full and tubes).
        """
        if ('incident_energy' not in info) \
                or (info['incident_energy'] == '350GeV'):
            e_lims = _e_lims_350
            tube_e_lims = _tube_e_lims_350
        else:
            e_lims = _e_lims_200
            tube_e_lims = _tube_e_lims_200

        return e_lims, tube_e_lims

    @staticmethod
    def _energy_label(e_dep, std=None):
        if std:
            return f"{int(round(e_dep))}$\\pm${int(round(std))} MeV"
        return f"{int(round(e_dep))} MeV"


class Event(Piece):
    """
    A single event.

        Attributes:
            hits: Data obtained from the hits file path.

            energy_vs_z: Energy vs. z histogram object, from calc.

            energy_vs_xy: See calc.

            numbers: See calc.

            run: The parent run that this event belongs to.
    """

    def __init__(self, hits_path, parent=None, out_dir=None, info=None):
        """

        :param hits_path: Path to particle hits data from Geant4.
        :type hits_path: str
        :param parent: The parent run that this event belongs to.
        :type parent: Run or None
        :param out_dir: Path to output directory.
        :type out_dir: str or None
        :param out_text: Path to output text. Defaults to
            out_dir/analysis.txt.
        :type out_text: str or None
        :param info: Specific parameters for this event.
        :type info: dict or None
        """
        super().__init__(hits_path, parent, out_dir, info)

        self.hits = None

        # For tagging Numbers entries.
        self.__tags = {'event': self.name}
        if self.parent:
            self.__tags['run'] = self.parent.name

        self.go()

    def _check_input_path(self, input_path):
        assert os.path.isfile(input_path), \
            "Event data path isn't a file."

    def _input_path2name(self, input_path):
        return _file2name(input_path)

    def go(self):
        """Read in data, make plots and do calculations."""
        self.hits = pandas.read_csv(self.input_path, skiprows=_skiprows)

        for histogram in self.hists.values():
            histogram.add_data(self.hits)
        self.numbers.add_data(self.hits, self.__tags)

        self.hists['energy_vs_z'].plot_single(
            energy_label=self._energy_label(self.numbers.get()['middle_e_dep'])
        )


class Run(Piece):
    """A collection of events with a given setup."""

    def __init__(self, events_path, out_dir=None, info=None):
        """

        :param events_path: Path to directory containing data from
            multiple events.
        :type events_path: str
        :param out_dir: Save things in here if this is given.
        :type out_dir: str or None
        :param info: Specific parameters for this run.
        :type info: dict or None
        """
        parent = None

        super().__init__(events_path, parent, out_dir, info)

        self.go()

    def go(self):
        """Create new analyzed Events and analyze this Run."""
        self._new_events()

        means, stds = self.numbers.append_mean_and_dev(
            mean_tags={'run': self.name, 'event': 'mean'},
            std_tags={'run': self.name, 'event': 'std_dev'}
        )
        self.numbers.save()

        energy_label = self._energy_label(
            means['middle_e_dep'], stds['middle_e_dep']
        )
        self.hists['energy_vs_z'].plot_means(energy_label=energy_label)
        self.hists['energy_vs_z'].plot_multi(energy_label=energy_label)

    def _new_events(self):
        """
        Create new analyzed Events and update this Run's analysis
        results.
        """
        for hits_path in self._event_paths():
            event = Event(hits_path, parent=self)

            for hist_name, hist in event.hists.items():
                self.hists[hist_name].add_results(hist.get())

            self.numbers.add_results(event.numbers.get())

    def _event_paths(self):
        """Get the paths to the event data files."""
        root, dirs, files = next(os.walk(self.input_path))
        for filename in files:
            if 'hits' in filename:
                yield os.path.join(root, filename)

    def _get_e_lims(self, info):
        # Update info here so that the Run name can be used.
        if 'incident_energy' not in info:
            # Parse the folder name.
            if '350gev' in self.name.lower():
                self.info['incident_energy'] = '350GeV'
            else:
                self.info['incident_energy'] = '200GeV'

        return super()._get_e_lims(info)

    def _check_input_path(self, input_path):
        assert os.path.isdir(input_path), \
            "Run directory path isn't actually a directory."

    def _input_path2name(self, input_path):
        return _dir2name(input_path)


def go(out_dir, *run_dirs):
    """
    Analyze a bunch of runs.

    :param out_dir: Where to save all the results, if at all.
    :type out_dir: str
    :param run_dirs: Locations of the run data directories.
    :type run_dirs: [str]
    :return:
    :rtype: None
    """
    if out_dir.lower() == 'none':
        out_dir = None

    for run_dir in run_dirs:
        Run(run_dir, out_dir, info=_run_params(run_dir))


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


def _weak_update(dict1, dict2):
    """
    Update dict1 with dict2 but with reversed priority. The existing
    keys and values in dict1 are kept.
    """
    new_dict1 = dict2.copy()
    new_dict1.update(dict1)
    dict1.update(new_dict1)
