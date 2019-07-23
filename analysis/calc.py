"""Energy calculations and histograms."""

import abc
import collections
import functools
import math
import os

import numpy
import pandas
import seaborn
from matplotlib import pyplot  #### plot  # TODO: Remove these plot __tags.

# Use seaborn's plotting styles.
seaborn.set()

_default_bin_density = 10
_default_dpi = 300
_default_length_units = 'mm'
_default_energy_units = 'MeV'

_default_z_lims = (-50, 80)
_default_tube_z_lims = (-35 / 2, 35 / 2)
_default_middle_z_lims = (-8 / 2, 8 / 2)

_default_xy_lims = 2 * ((-20 / 2, 20 / 2),)
_default_xy_hist_z_lims = _default_middle_z_lims

_image_format = 'jpg'
_data_format = 'csv'
_linewidth = 0.5


class Calc(abc.ABC):
    """
    A base for calculations and presentations of results a certain
    kind.

        Attributes:
            resultss: A collection of the results. It is a pandas DataFrame
            in this base class, but that can change.

            piece: Piece of analysis that this calculation is for.

            save_dir: Save output in this folder.

            save_name: A name for the saved files.
    """

    def __init__(self, piece):
        self._piece = piece
        self.resultss = pandas.DataFrame()

        if self._piece.out_dir:
            # Save when Piece has an output directory.
            self.save_dir = self._piece.out_dir
            self.save_name = self._piece.name
        else:
            self.save_dir = None
            self.save_name = None

    def get(self, i=None):
        """
        Get the ith set of results, or if i isn't given, make sure
        resultss contains only one entry and return it.

        :param i: The index of the results to get.
        :type i: int or None
        :return: Results.
        :rtype:
        """
        if i is None:
            self._assert_single_entry()

            return self.resultss.loc[0]

        return self.resultss.loc[i]

    def add_data(self, data):
        """
        Calculate results from data and keep them.

        :param data: Hits data.
        :type data: pandas.DataFrame
        :return:
        :rtype:
        """
        self.add_results(self._data2results(data))

    def add_results(self, results):
        """
        Keep some already calculated results (after making sure they are
        good results).

        :param results: Results to add.
        :type results: pandas.Series, dict, collections.OrderedDict
        :return:
        :rtype: None
        """
        self._check_results(results)
        self.resultss = self.resultss.append(results, ignore_index=True)

    @abc.abstractmethod
    def _data2results(self, data):
        """
        Calculate results from hits data and return them.

        :param data: Hits data.
        :return: Calc results.
        """

    @abc.abstractmethod
    def _check_results(self, results):
        """Make sure some calculated results are good to keep."""

    def _assert_single_entry(self):
        """Make sure only one set of results has been added."""
        assert len(self.resultss) == 1, \
            "This Calc container has more than one set of results."

    def _assert_multiple_entries(self):
        """Make sure more than one set of results has been added."""
        assert len(self.resultss) > 1, \
            "This Calc container does not have more than one set of" \
            " results."


class Numbers(Calc):
    """
    The numbers don't lie.

    This class manages calculations of different numbers from the data
    and writes them to files.
    """

    def add_data(self, data, tags=None):
        self.add_results(self._data2results(data, tags))

    def append_mean_and_dev(self, mean_tags=None, std_tags=None):
        """
        Append mean and standard deviation entries to the collection of
        results.

        :param mean_tags: Tags to use for the new means entry.
        :type mean_tags: dict or collections.OrderedDict or None
        :param std_tags: Tags to use for the new standard dev entry.
        :type std_tags: dict or collections.OrderedDict or None
        :return:
        :rtype:
        """
        # Pandas' std function gives NANs for DataFrames with only one
        # entry.
        self._assert_multiple_entries()

        means = self.resultss.mean()
        stds = self.resultss.std()

        if mean_tags:
            means = means.append(pandas.Series(mean_tags))
        if std_tags:
            stds = stds.append(pandas.Series(std_tags))

        self.add_results(means)
        self.add_results(stds)

    def save(self):
        """
        Save results to a csv file if save_dir is set.

        :return: Path to saved results file.
        :rtype: str
        """
        if not self.save_dir:
            return None

        filename = '.'.join((self.save_name, _data_format))
        filepath = os.path.join(self.save_dir, filename)
        self.resultss.to_csv(filepath, index=False)

        return filepath

    def _data2results(self, data, tags=None):
        """

        :param data:
        :type data:
        :param tags: Extra result values for tagging an entry (e.g. with
            its event and run).
        :type tags: dict or collections.OrderedDict or None
        :return:
        :rtype:
        """
        new_results = collections.OrderedDict()

        if tags:
            new_results.update(tags)

        new_results.update(self.__split_z(data))
        new_results.update(self.__tubes(data))

        # TODO: Remove these.

        # new_results.update(self._piece.__tags)

        # self.resultss.append(new_results, ignore_index=True)
        # self._append(new_results)

        return new_results

    def _check_results(self, results):
        assert len(results) <= 9, "Tried to add funny results."

    @staticmethod
    def __split_z(data):
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

    @staticmethod
    def __tubes(data):
        """TODO: Document and organize this."""
        offset = 7.5 / 4
        data = data[numpy.logical_and(_default_middle_z_lims[0] < data.z, data.z < _default_middle_z_lims[1])]

        right = data.z > 0
        left = data.z <= 0
        top_right = numpy.logical_and(right, data.y > offset)
        bottom_right = numpy.logical_and(right, data.y < offset)
        top_left = numpy.logical_and(left, data.y > -offset)
        bottom_left = numpy.logical_and(left, data.y < -offset)
        top_right_sum = data.energy_deposit[top_right].sum()
        top_left_sum = data.energy_deposit[top_left].sum()
        bottom_left_sum = data.energy_deposit[bottom_left].sum()
        bottom_right_sum = data.energy_deposit[bottom_right].sum()
        return {'top_right': top_right_sum, 'top_left': top_left_sum, 'bottom_left': bottom_left_sum, 'bottom_right': bottom_right_sum}


class Histogram(Calc):
    """
    Base histogram.

    TODO: Close figures.

        Attributes:
            bin_density: E.g. bins per mm.

            dpi: Dots per inch for images.

            title: The plot title.
    """

    def __init__(
            self,
            piece,
            title=None,
            bin_density=None,
            dpi=None,
            energy_units=None,
            length_units=None
    ):
        super().__init__(piece)

        # Attributes with defaults.
        self.title = self._default_title or title
        self.bin_density = bin_density or _default_bin_density
        self.dpi = dpi or _default_dpi
        self.energy_units = energy_units or _default_energy_units
        self.length_units = length_units or _default_length_units

    def save_fig(self, fig, file_suffix):
        """
        Save a figure to a new file if save_dir is set.

        :param fig: The figure to save.
        :type fig: figure
        :param file_suffix: Added to the end of the saved filename.
        :type file_suffix: str or None
        :return: File path of the saved figure.
        :rtype: str
        """
        if not self.save_dir:
            return None

        if file_suffix:
            save_name = '-'.join((self.save_name, file_suffix))
        else:
            save_name = self.save_name

        filename = '.'.join((save_name, _image_format))
        filepath = os.path.join(self.save_dir, filename)

        fig.savefig(
            filepath, dpi=self.dpi, format=_image_format, bbox_inches='tight'
        )

        return filepath

    def _to_density(self, sums):
        """
        Convert values from units of energy/bin to units of
        energy/distance.

        :param sums: Values to convert.
        :type sums: numpy.ndarray
        :return: Converted values.
        :rtype: numpy.ndarray
        """
        return sums * self.bin_density

    @property
    def _default_title(self):
        """Default plot title."""
        raise NotImplementedError

    @abc.abstractmethod
    def plot_single(self, i=None, save=True):
        """
        Plot the ith set of sums. If i isn't given, assume there is only
        set of sums and plot it.

        :param save: Save to file with the standard filename if True.
        :type save: bool
        :param i: Index of the sums.
        :type i: int
        :return: New figure and axis/axes.
        :rtype: tuple
        """

    @abc.abstractmethod
    def plot_means(self):
        """
        Plot the mean, and maybe the standard deviation, of all results.

        :return: New figure and axis/axes.
        :rtype: tuple
        """

    @abc.abstractmethod
    def _new_fig_and_axes(self):
        """
        Make a new labeled figure and axis/axes.

        :return: The figure and axis/axes.
        :rtype: tuple
        """

    def plot_single_sums(self, sums):
        """
        Plot a lone set of histogram sums.
        TODO: Remove.

        :param sums: Histogram sums.
        :type sums: numpy.ndarray
        :return: New figure and axis/axes.
        :rtype: tuple
        """


class EnergyVsZ(Histogram):
    """
    Histograms of energy deposit vs. z.

        Attributes:
            e_lims: Overall energy/y-axis limits.

            tube_e_lims: Tube energy/y-axis limits.

            z_lims: Overall z limits.

            tube_z_lims: Tube z limits.
    """

    _default_title = 'Energy vs. z.'

    def __init__(
            self,
            piece,
            e_lims,
            tube_e_lims,
            title=None,
            z_lims=None,
            tube_z_lims=None
    ):
        super().__init__(piece, title)

        self.__e_lims = e_lims
        self.__tube_e_lims = tube_e_lims

        # Attributes with defaults.
        self.__z_lims = z_lims or _default_z_lims
        self.__tube_z_lims = tube_z_lims or _default_tube_z_lims

        self.__bins = _make_bins(
            self.__z_lims[0], self.__z_lims[1], self.bin_density
        )
        self.__bin_mids = _make_bin_midpoints(self.__bins)

    def plot_single(self, i=None, save=True):
        fig, ax, ax_middle = self._new_fig_and_axes()
        fig.suptitle(self.title)

        plate_kwargs = {'linewidth': _linewidth, 'color': 'blue'}
        tube_kwargs = {'linewidth': _linewidth, 'color': 'purple'}
        _split_plot(
            x=self.__bin_mids,
            y=self._to_density(self.get(i)),
            split_lims=self.__tube_z_lims,
            ax1=ax_middle,
            ax2=ax,
            kwargs1=plate_kwargs,
            kwargs2=tube_kwargs
        )

        if save:
            self.save_fig(fig, file_suffix=None)

    def plot_means(self):
        fig, ax, ax_middle = self._new_fig_and_axes()
        fig.suptitle(' '.join((self.title, '(Mean)')))

        plate_kwargs = {'linewidth': _linewidth, 'color': 'blue'}
        tube_kwargs = {'linewidth': _linewidth, 'color': 'purple'}
        _split_plot(
            x=self.__bin_mids,
            y=self._to_density(numpy.mean(self.resultss, axis=0)),
            split_lims=self.__tube_z_lims,
            ax1=ax_middle,
            ax2=ax,
            kwargs1=plate_kwargs,
            kwargs2=tube_kwargs
        )

        self.save_fig(fig, file_suffix=None)

    def _data2results(self, data):
        return pandas.Series(
            numpy.histogram(
                data.z, bins=self.__bins, weights=data.energy_deposit
            )[0],
            index=self.__bin_mids
        )

    def _check_results(self, sums):
        assert len(sums) == (len(self.__bins) - 1), \
            "Tried to add funny histogram sums."

    def _new_fig_and_axes(self):
        fig, ax = pyplot.subplots()
        ax_middle = ax.twinx()

        ax.set_xlabel(f'z ({self.length_units})')
        ax.set_ylabel(
            'E Dep Density in Plate Cals'
            f' ({self.energy_units} / {self.length_units})'
        )
        ax_middle.set_ylabel(
            'E Dep Density in Tube Cals'
            f' ({self.energy_units} / {self.length_units})'
        )

        for axes, e_limits in zip(
                (ax, ax_middle), (self.__e_lims, self.__tube_e_lims)
        ):
            axes.set_xlim(self.__z_lims)
            axes.set_ylim(e_limits)

        return fig, ax, ax_middle


class EnergyVsXY(Histogram):
    """Histogram of energy deposit vs. x and y."""

    _default_title = 'Energy vs. x and y.'

    def __init__(self, piece, title=None, xy_lims=None, z_lims=None):
        super().__init__(piece, title)

        # Attributes with defaults.
        self.__xy_lims = xy_lims or _default_xy_lims
        self.__z_lims = z_lims or _default_xy_hist_z_lims

        self.__binss = tuple(
            _make_bins(limits[0], limits[1], self.bin_density)
            for limits in self.__xy_lims
        )
        self.__bin_midss = tuple(
            _make_bin_midpoints(bins) for bins in self.__binss
        )

    def plot_single(self, i=None, save=True):
        fig, ax = self._new_fig_and_axes()
        fig.suptitle(self.title)

        grid_x, grid_y = numpy.meshgrid(*self.__binss)
        ax.pcolormesh(grid_x, grid_y, self.get(i))

        self.save_fig(fig, file_suffix=None)

    def plot_means(self):
        pass

    def _data2results(self, data):
        sliced = data[numpy.logical_and(
            self.__z_lims[0] <= data.z, data.z <= self.__z_lims[1]
        )]
        return numpy.histogram2d(
            sliced.x, sliced.y, bins=self.__binss, weights=sliced.energy_deposit
        )[0]

    def _check_results(self, sums):
        """Check that `sums` is the correct shape for the bins."""
        assert numpy.all(numpy.equal(
            sums.shape,
            tuple(len(bins) - 1 for bins in self.__binss)
        ))

    def _new_fig_and_axes(self):
        fig, ax = pyplot.subplots()

        ax.set_xlabel(f'x ({self.length_units})')
        ax.set_ylabel(f'y ({self.length_units})')

        ax.set_xlim(self.__xy_lims[0])
        ax.set_ylim(self.__xy_lims[1])

        return fig, ax


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
    TODO: Remove this.

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
