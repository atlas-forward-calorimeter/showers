"""Energy calculations and histograms."""

# TODO: Consider moving save paths to the save methods;
#  but this isn't necessary.

import abc
import collections
import functools
import math
import os

import numpy
import pandas
from matplotlib import pyplot

_default_bin_density = 10
_default_length_units = 'mm'
_default_energy_units = 'MeV'
_default_dpi = 300

_default_z_lims = (-50, 80)
_default_tube_z_lims = (-35 / 2, 35 / 2)
_default_middle_z_lims = (-8 / 2, 8 / 2)

_default_xy_lims = 2 * ((-20 / 2, 20 / 2),)
_default_xy_hist_z_lims = _default_middle_z_lims

_image_format = 'jpg'
_data_format = 'csv'
_linewidth = 0.5

SplitZResults = collections.namedtuple('SplitZResults', 'full tubes middle')


class Calc(abc.ABC):
    """
    A base for calculations and presentations of results a certain
    kind.

    Attributes:
        resultss: A collection of the results.
    """

    def __init__(self):
        self.resultss = []

    def get(self, i=None):
        """
        Get the ith set of results kept in resultss. By default, make
        sure that resultss contains only one entry, then return it.

        :param i: The index of the results to get.
        :type i: int or None
        :return: Results.
        :rtype:
        """
        if i is None:
            self._assert_single_entry()

            return self.resultss[0]

        return self.resultss[i]

    def _assert_single_entry(self):
        """Make sure only one set of results has been added."""
        assert len(self.resultss) == 1, \
            "This Calc container has more than one set of results."

    @abc.abstractmethod
    def add_data(self, data):
        """
        Calculate results from data and keep them.

        :param data: Hits data.
        :type data: pandas.DataFrame
        :return: Histogram sums.
        :rtype: numpy.ndarray
        """

    @abc.abstractmethod
    def add_results(self, results):
        """
        Add some already calculated results (after optionally making
        sure they are good results).

        :param results: Results to add.
        :type results: numpy.ndarray
        :return:
        :rtype: None
        """


class Histogram(Calc):
    """
    Base histogram class.

    Attributes:
        resultss: Collection of all histogram sums.
    """
    # TODO: Close figures.

    def __init__(self, save_name=None, save_dir=None, bin_density=None,
                 dpi=None, energy_units=None, length_units=None):
        """

        :param bin_density: Bins per length.
        :type bin_density: float
        :param energy_units: Energy dimensions.
        :type energy_units: str or None
        :param length_units: Length dimensions.
        :type length_units: str or None
        :param save_name: Name of saved figure files. Suffixes are added.
            Figures are only saved when save_dir is given.
        :type save_name: str or None
        :param save_dir: Save figures to this directory, and only save
            when this is given (since Event and Run also work that way).
        :type save_dir: str or None
        :param dpi: Pixels per inch for saved figures.
        :type dpi: int or None
        """
        super().__init__()

        # defaults
        self.bin_density = _default_bin_density if bin_density is None \
            else bin_density
        self.dpi = _default_dpi if dpi is None else dpi
        self.energy_units = _default_energy_units if energy_units is None \
            else energy_units
        self.length_units = _default_length_units if length_units is None \
            else length_units

        self.save_name, self.save_dir = _default_save_params(
            save_name, save_dir
        )
        ####

    def save_fig(self, fig, file_suffix):
        """
        Save a figure to a new file if save_name is defined.

        :param fig: The figure to save.
        :type fig: figure
        :param file_suffix: Added to the end of the saved filename.
        :type file_suffix: str or None
        :return: File path of the saved figure.
        :rtype: str
        """
        # Only save when save_dir is given.
        if self.save_dir is None:
            return None

        if file_suffix is None:
            save_name = self.save_name
        else:
            save_name = '_'.join((self.save_name, file_suffix))

        filename = '.'.join((save_name, _image_format))
        filepath = os.path.join(self.save_dir, filename)

        fig.savefig(
            filepath, dpi=self.dpi, format=_image_format, bbox_inches='tight'
        )

        return filepath

    def to_density(self, array):
        """
        Convert values from units of energy/bin to units of
        energy/distance.

        :param array: Values to convert.
        :type array: numpy.ndarray
        :return: Converted values.
        :rtype: numpy.ndarray
        """
        return array * self.bin_density

    def plot_single(self, i=None):
        """
        Plot the ith set of sums kept in sumss. By default, make sure
        that sumss has only one entry, then plot it.

        :param i: The index.
        :type i: int
        :return: New figure and axis/axes.
        :rtype: tuple
        """
        if i is None:
            self._assert_single_entry()

            return self.plot_single_sums(self.resultss[0])

        return self.plot_single_sums(self.resultss[i])

    @abc.abstractmethod
    def plot_single_sums(self, sums):
        """
        Plot a lone set of histogram sums.

        :param sums: Histogram sums.
        :type sums: numpy.ndarray
        :return: New figure and axis/axes.
        :rtype: tuple
        """

    @abc.abstractmethod
    def new_fig_and_axes(self):
        """
        Make a new labeled figure and axis/axes.

        :return: The figure and axis/axes.
        :rtype: tuple
        """


class EnergyVsZ(Histogram):
    """Histograms of energy deposit vs. z."""

    def __init__(self, e_lims, tube_e_lims, e_units=None, z_lims=None,
                 tube_z_lims=None, z_units=None, bin_density=None, title=None,
                 save_name=None, save_dir=None, dpi=None):
        """

        :param save_name: Filename prefix of saved figures.
        :type save_name: str
        :param save_dir: Directory of saved figures.
        :type save_dir: str
        :param dpi: Pixels per inch for saved figure.
        :type dpi: int or None
        :param z_lims: (lower z limit, upper z limit)
        :type z_lims: (float, float) or None
        :param tube_z_lims: (lower z limit, upper z limit)
        :type tube_z_lims: (float, float) or None
        :param z_units: z dimensions label.
        :type z_units: str or None
        :param e_lims: (lower energy limit, upper energy limit)
        :type e_lims: (float, float)
        :param tube_e_lims: (lower energy limit, upper energy limit)
        :type tube_e_lims: (float, float)
        :param e_units: Energy dimensions label.
        :type e_units: str or None
        :param bin_density: Bins per z.
        :type bin_density: float or None
        :param title: Plot title. Suffixes are added.
        :type title: str or None
        """
        super().__init__(
            save_name, save_dir, bin_density, dpi, e_units, z_units
        )

        self._e_lims = e_lims
        self._tube_e_lims = tube_e_lims

        # defaults
        self._z_lims = _default_z_lims if z_lims is None else z_lims
        self._tube_z_lims = _default_tube_z_lims if tube_z_lims is None \
            else tube_z_lims

        self._title = 'Energy vs. z.' if title is None else title
        ####

        self._bins = _make_bins(
            self._z_lims[0], self._z_lims[1], self.bin_density
        )

    def add_data(self, data):
        sums, _ = numpy.histogram(
            data.z, bins=self._bins, weights=data.energy_deposit
        )
        self.resultss.append(sums)

        return sums

    def add_results(self, sums):
        assert len(sums) == (len(self._bins) - 1), \
            "Tried to add funny sums."

        self.resultss.append(sums)

    def plot_single_sums(self, sums):
        bin_mids = _make_bin_midpoints(self._bins)

        fig, ax, ax_middle = self.new_fig_and_axes()
        fig.suptitle(self._title)

        plot_kwargs = {'linewidth': _linewidth, 'color': 'purple'}
        _split_plot(
            bin_mids, self.to_density(sums), self._tube_z_lims,
            ax_middle, ax, plot_kwargs, plot_kwargs
        )

        self.save_fig(fig, file_suffix=None)

    def new_fig_and_axes(self):
        fig, ax = pyplot.subplots()
        ax_middle = ax.twinx()

        ax.set_xlabel(f'z ({self.length_units})')
        ax.set_ylabel(
            'Energy Deposit Density'
            f' ({self.energy_units} / {self.length_units})'
        )

        for axes, e_limits in zip(
                (ax, ax_middle), (self._e_lims, self._tube_e_lims)
        ):
            axes.set_xlim(self._z_lims)
            axes.set_ylim(e_limits)

        return fig, ax, ax_middle


class EnergyVsXY(Histogram):
    """Histogram of energy deposit vs. x and y."""

    def __init__(self, title=None, xy_lims=None, z_lims=None, save_name=None,
                 save_dir=None, bin_density=None, dpi=None, energy_units=None,
                 length_units=None):
        """

        :param title: Plot title. Suffixes are added.
        :type title: str
        :param xy_lims: ((x lower limit, x upper limit),
            (y lower limit, y upper limit))
        :type xy_lims: ((float, float), (float, float))
        :param z_lims: (lower limit, upper limit)
        :type z_lims: (float, float)
        :param save_name:
        :type save_name:
        :param save_dir:
        :type save_dir:
        :param bin_density:
        :type bin_density:
        :param dpi:
        :type dpi:
        :param energy_units:
        :type energy_units:
        :param length_units:
        :type length_units:
        """
        super().__init__(save_name, save_dir, bin_density, dpi, energy_units,
                         length_units)

        # defaults
        self._xy_lims = _default_xy_lims if xy_lims is None else xy_lims
        self._z_lims = _default_xy_hist_z_lims if z_lims is None else z_lims

        self._title = 'Energy vs. x & y.' if title is None else title
        #

        self._binss = tuple(
            _make_bins(limits[0], limits[1], self.bin_density)
            for limits in self._xy_lims
        )

    def add_data(self, data):
        sliced = data[numpy.logical_and(
            self._z_lims[0] <= data.z, data.z <= self._z_lims[1]
        )]
        sums, _, _ = numpy.histogram2d(
            sliced.x, sliced.y, bins=self._binss, weights=sliced.energy_deposit
        )
        self.resultss.append(sums)

        return sums

    def add_results(self, sums):
        if self.resultss:
            assert numpy.all(numpy.equal(
                sums.shape,
                tuple(len(bins) - 1 for bins in self._binss)
            ))

        self.resultss.append(sums)

    def plot_single_sums(self, sums):
        # bin_midss = (_make_bin_midpoints(bins) for bins in self._binss)

        fig, ax = self.new_fig_and_axes()
        fig.suptitle(self._title)

        grid_x, grid_y = numpy.meshgrid(*self._binss)
        ax.pcolormesh(grid_x, grid_y, sums)

        self.save_fig(fig, file_suffix=None)

    def new_fig_and_axes(self):
        fig, ax = pyplot.subplots()

        ax.set_xlabel(f'x ({self.length_units})')
        ax.set_ylabel(f'y ({self.length_units})')

        ax.set_xlim(self._xy_lims[0])
        ax.set_ylim(self._xy_lims[1])

        return fig, ax


class Numbers(Calc):
    """
    The numbers don't lie.

    This class manages calculations of different numbers from the data
    and writes them to files.

    Attributes:
        resultss: Collection of all results. Results are dicts.
    """

    def __init__(self, save_name=None, save_dir=None):
        """

        :param save_name: See Histogram.
        :type save_name:
        :param save_dir: See Histogram.
        :type save_dir:
        """
        super().__init__()

        self.dataframe_resultss = None

        self._save_name, self._save_dir = _default_save_params(
            save_name, save_dir
        )

    def add_data(self, data, tags=None):
        """

        :param data:
        :type data:
        :param tags: Extra result values for tagging a result entry.
        :type tags: dict or None
        :return:
        :rtype:
        """
        new_results = collections.OrderedDict()

        if tags is not None:
            new_results.update(tags)

        new_results.update(self._split_z(data))
        self.resultss.append(new_results)

        return new_results

    def add_results(self, results):
        assert len(results) <= 5, "Tried to add funny results."

        self.resultss.append(results)

    def append_mean_and_dev(self, mean_tags=None, dev_tags=None):
        """
        Append mean and standard deviation entries to the collection of
        results.

        :param mean_tags: Tags to use for the new means entry.
        :type mean_tags: dict
        :param dev_tags: Tags to use for the new dev entry.
        :type dev_tags: dict
        :return:
        :rtype:
        """
        assert self.resultss, \
            "Can't take the mean and standard deviation of" \
            " nonexistent `resultss`."

        df_resultss = pandas.DataFrame(self.resultss)

        means = df_resultss.mean().to_dict()
        devs = df_resultss.std().to_dict()

        if mean_tags:
            means.update(mean_tags)
        if dev_tags:
            devs.update(dev_tags)

        self.add_results(means)
        self.add_results(devs)

    def save(self):
        """
        Save results to a file that can be read by pandas, Excel, etc.
        if _save_dir is set.

        :return: Path to saved results file.
        :rtype: str
        """
        if self._save_dir is None:
            return None

        filename = '.'.join((self._save_name, _data_format))
        filepath = os.path.join(self._save_dir, filename)
        pandas.DataFrame(self.resultss).to_csv(filepath, index=False)

        return filepath

    @staticmethod
    def _split_z(data):
        """
        Calculate energy deposit mean and std dev on sections of
        data split up along z.

        :param data: Hits data.
        :type data: pandas.DataFrame
        :return: Calculation results.
        :rtype: dict
        """
        tube_indices = numpy.logical_and(
            _default_tube_z_lims[0] <= data.z,
            data.z <= _default_tube_z_lims[1]
        )
        middle_indices = numpy.logical_and(
            _default_middle_z_lims[0] <= data.z,
            data.z <= _default_middle_z_lims[1]
        )
        return {
            'full_e_dep':
                data.energy_deposit.sum(),
            'tubes_e_dep':
                data.energy_deposit[tube_indices].sum(),
            'middle_e_dep':
                data.energy_deposit[middle_indices].sum(),
        }


@functools.lru_cache()
def _make_bins(start, end, bin_density):
    """
    Make histogram bins with a given bin density.

    Since the bin density is fixed, the bins start exactly at `start`
    and end at a rounded up value near `end`.

    :param start: Exact start of the first bin.
    :type start: float
    :param end: Lower limit of the end of last bin.
    :type end: float
    :param bin_density: Exact bin density.
    :type bin_density: float
    :return: Bin endpoints.
    :rtype: numpy.ndarray
    """
    assert (end > start)

    # Round the end limit so the bin size is correct.
    num_bins = math.ceil(bin_density * (end - start))
    new_end = start + num_bins / bin_density

    return numpy.linspace(start, new_end, num_bins + 1)


# @functools.lru_cache()
def _make_bin_midpoints(bins):
    """Get midpoints of histogram bins."""
    return (bins[:-1] + bins[1:]) / 2


def _range2indices(array, range_lims):
    """
    Find the values of an array that are inclusively within range_lims.

    :param array: Array of values that are compared to range_lims.
    :type array: numpy.ndarray
    :param range_lims: (lower limit, upper limit)
    :type range_lims: (float, float)
    :return: Boolean array locating the indices of the values within the
        range.
    :rtype: numpy.ndarray
    """
    return numpy.logical_and(range_lims[0] <= array, array <= range_lims[1])


def _default_save_params(save_name, save_dir):
    """
    Require save_name if save_dir is given.

    :param save_name: Name of file to save.
    :type save_name: str or None
    :param save_dir: Directory to save to.
    :type save_dir: str or None
    :return: save_name, save_dir
    :rtype: (str, str)
    """
    if save_dir is not None:
        assert save_name is not None, (
            "Please give a save filename if you would like to"
            " save figures."
        )

    return save_name, save_dir


def _split_plot(x, y, split_lims, ax1, ax2, kwargs1=None, kwargs2=None):
    """
    Split up data and plot it on two axes.

    The data with x values inside split_lims (inclusive) goes on ax1,
    and the rest goes on ax2.

    :param x: x values.
    :type x: 1d array
    :param y: y values.
    :type y: 1d array
    :param split_lims: Split data based on x values within or outside of
        these limits.
    :type split_lims: (float, float)
    :param ax1: Data inside split_lims (inclusive) goes on this.
    :type ax1: axes
    :param ax2: Data outside split_lims goes on this.
    :type ax2: axes
    :param kwargs1: Keyword arguments passed to the ax1 plot.
    :type kwargs1: dict or None
    :param kwargs2: Keyword arguments passed to the ax2 plot.
    :type kwargs2: dict or None
    :return:
    :rtype:
    """
    if kwargs1 is None:
        kwargs1 = {}
    if kwargs2 is None:
        kwargs2 = {}

    inside_indices = _range2indices(x, split_lims)
    outside_indices = numpy.logical_not(inside_indices)

    ax1.plot(x[inside_indices], y[inside_indices], **kwargs1)
    ax2.plot(x[outside_indices], y[outside_indices], **kwargs2)
