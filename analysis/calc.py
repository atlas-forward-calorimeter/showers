"""Energy calculations and histograms."""

import abc
import collections
import functools
import math
import os

import numpy
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

_format = 'jpg'
_linewidth = 0.5

SplitZResults = collections.namedtuple('SplitZResults', 'full tubes middle')


class Histogram(abc.ABC):
    """
    Base histogram class.

    Attributes:
        sumss: Collection of all histogram sums.
    """

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
        :type save_name: str or None
        :param save_dir: Save figures to this directory.
        :type save_dir: str or None
        :param dpi: Pixels per inch for saved figures.
        :type dpi: int or None
        """
        self.sumss = []

        # defaults

        self.bin_density = _default_bin_density if bin_density is None \
            else bin_density
        self.dpi = _default_dpi if dpi is None else dpi
        self.energy_units = _default_energy_units if energy_units is None \
            else energy_units
        self.length_units = _default_length_units if length_units is None \
            else length_units

        if save_name is None:
            assert save_dir is None, (
                "Please give a save filename if you would like to"
                " save figures."
            )
            self.save_name = None
            self.save_dir = None
        else:
            self.save_name = save_name
            if save_dir is None:
                self.save_dir = ''
            else:
                self.save_dir = save_dir

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
        if self.save_name is None:
            return None

        if file_suffix is None:
            save_name = self.save_name
        else:
            save_name = '_'.join((self.save_name, file_suffix))

        filename = '.'.join((save_name, _format))
        filepath = os.path.join(self.save_dir, filename)

        fig.savefig(
            filepath, dpi=self.dpi, format=_format, bbox_inches='tight'
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

    @abc.abstractmethod
    def add_data(self, data):
        """
        Calculate sums from data and keep them.

        :param data: Hits data.
        :type data: pandas.DataFrame
        :return: Histogram sums.
        :rtype: numpy.ndarray
        """

    @abc.abstractmethod
    def add_sums(self, sums):
        """
        Add some already calculated sums after making sure they are good
        sums.

        :param sums: Histogram sums to add.
        :type sums: numpy.ndarray
        :return:
        :rtype: None
        """

    @abc.abstractmethod
    def plot_single(self, i=0):
        """
        Plot the ith set of kept sums.

        :param i: Index, defaults to 0.
        :type i: int
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

        self._bins = make_bins(
            self._z_lims[0], self._z_lims[1], self.bin_density
        )

    def add_data(self, data):
        sums, _ = numpy.histogram(
            data.z, bins=self._bins, weights=data.energy_deposit
        )
        self.sumss.append(sums)

        return sums

    def add_sums(self, sums):
        assert len(sums) == (len(self._bins) - 1), \
            "Tried to add funny sums."

        self.sumss.append(sums)

    def plot_single(self, i=0):
        sums = self.sumss[i]

        bin_mids = make_bin_midpoints(self._bins)

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
            make_bins(limits[0], limits[1], self.bin_density)
            for limits in self._xy_lims
        )

    def add_data(self, data):
        sliced = data[numpy.logical_and(
            self._z_lims[0] <= data.z, data.z <= self._z_lims[1]
        )]
        sums, _, _ = numpy.histogram2d(
            sliced.x, sliced.y, bins=self._binss, weights=sliced.energy_deposit
        )
        self.sumss.append(sums)

        return sums

    def add_sums(self, sums):
        assert numpy.all(numpy.equal(
            sums.shape,
            (len(self._binss[0]), len(self._binss[1]))
        ))

        self.sumss.append(sums)

    def plot_single(self, i=0):
        sums = self.sumss[i]

        bin_midss = (make_bin_midpoints(bins) for bins in self._binss)

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


class SplitZ:
    """
    Calculation of energy deposit mean and std dev on sections of
    data split up along z.

    Attributes:
        resultss: Collection of all results.
    """

    def __init__(self):
        self.resultss = []

    def add_data(self, data):
        """
        Calculate results from energy deposit data and keep them.

        :param data: Hits data.
        :type data: pandas.DataFrame
        :return: New results. (total, tubes, middle tubes)
        :rtype: (float, float, float)
        """
        tube_indices = numpy.logical_and(
            _default_tube_z_lims[0] <= data.z,
            data.z <= _default_tube_z_lims[1]
        )
        middle_indices = numpy.logical_and(
            _default_middle_z_lims[0] <= data.z,
            data.z <= _default_middle_z_lims[1]
        )
        results = SplitZResults(*tuple(
            numpy.sum(data.energy_deposit[indices])
            for indices in (slice(None), tube_indices, middle_indices)
        ))
        self.resultss.append(results)

        return results

    def add_results(self, results):
        """
        Add some already calculated results.

        :param results: Results to add. See add_data.
        :type results: tuple
        :return:
        :rtype:
        """
        assert len(results) == 3, "Tried to add funny results."

        self.resultss.append(results)

    def get(self, i=0):
        """Get the ith set of results."""
        return self.resultss[i]


@functools.lru_cache()
def make_bins(start, end, bin_density):
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
def make_bin_midpoints(bins):
    """Get midpoints of histogram bins."""
    return (bins[:-1] + bins[1:]) / 2


def range2indices(array, range_lims):
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

    inside_indices = range2indices(x, split_lims)
    outside_indices = numpy.logical_not(inside_indices)

    ax1.plot(x[inside_indices], y[inside_indices], **kwargs1)
    ax2.plot(x[outside_indices], y[outside_indices], **kwargs2)
