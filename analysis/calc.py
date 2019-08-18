"""Energy calculations and histograms."""

import abc
import collections
import datetime
import functools
import itertools
import math
import os

import numpy
import pandas
import seaborn
from matplotlib import pyplot

# Use seaborn's plotting styles.
seaborn.set()

# Energy-z histogram limits for the 8 by 4 plate arrangements. The
# limits for the 7 by 7 and 6 by 6 arrangements are shifted towards
# -z in increments of length equal to the `_plate_separation`.
z_lims = {(8, 4): (-65, 95)}
_plate_separation = 4 + 3.5
for i, plates in enumerate(((7, 5), (6, 6)), start=1):
    z_lims[plates] = tuple(
        # Subtract multiples of the `_plate_separation` from the 8 by 4
        # limits.
        lim - _plate_separation * i for lim in z_lims[(8, 4)]
    )

_default_bin_density = 10
_default_dpi = 300
_default_length_units = 'mm'
_default_energy_units = 'MeV'

_default_tube_z_lims = (-35 / 2, 35 / 2)

# Convert the _default_middle_z_lims to integers so that the plot
# tick labels are displayed without decimals.
_default_middle_z_lims = tuple(int(lim) for lim in (-8 / 2, 8 / 2))

_default_xy_lims = 2 * ((-3, 3),)
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

    def __repr__(self):
        return repr(self.resultss)

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
        :type results: pandas.Series or dict or collections.OrderedDict
        :return:
        :rtype: None
        """
        self._check_results(results)
        self.resultss = _ordered_append(self.resultss, results)

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

        return save_dataframe(self.resultss, filepath)

    def _assert_single_entry(self):
        """Make sure only one set of results has been added."""
        assert len(self.resultss) == 1, \
            "This Calc container does not have a single entry."

    def _assert_multiple_entries(self):
        """Make sure more than one set of results has been added."""
        assert len(self.resultss) > 1, \
            "This Calc container does not have more than one set of" \
            " results."

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
            means = pandas.Series(mean_tags).combine_first(means)
        if std_tags:
            stds = pandas.Series(std_tags).combine_first(stds)

        self.add_results(means)
        self.add_results(stds)

        return means, stds

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
            'tube_e_dep':
                data.energy_deposit[tube_indices].sum(),
            'middle_e_dep':
                data.energy_deposit[middle_indices].sum(),
        }

    @staticmethod
    def __tubes(data):
        """TODO: Document and organize this."""
        offset = 7.5 / 4
        # TODO: Change default to self (probably don't need to do this).
        data = data[numpy.logical_and(
            _default_middle_z_lims[0] < data.z,
            data.z < _default_middle_z_lims[1]
        )]

        bottom = data.y <= 0
        top = data.y > 0

        bottom_left = numpy.logical_and(bottom, data.x <= -offset)
        bottom_right = numpy.logical_and(bottom, data.x > -offset)
        top_left = numpy.logical_and(top, data.x <= offset)
        top_right = numpy.logical_and(top, data.x > offset)

        bottom_left_sum = data.energy_deposit[bottom_left].sum()
        bottom_right_sum = data.energy_deposit[bottom_right].sum()
        top_left_sum = data.energy_deposit[top_left].sum()
        top_right_sum = data.energy_deposit[top_right].sum()

        return {
            'top_right': top_right_sum,
            'top_left': top_left_sum,
            'bottom_left': bottom_left_sum,
            'bottom_right': bottom_right_sum
        }


class Histogram(Calc):
    """
    Base histogram.

        Attributes:
            bin_density: E.g. bins per mm.

            dpi: Dots per inch for images.

            title: The plot title.

    TODO: Add the option to close all the plot figures created.
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
        self.title = title or self._default_title
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
    def _make_fig_and_axes(self):
        """
        Make a new labeled figure and axis/axes.

        :return: The figure and axis/axes.
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

            middle_z_lims: Tube middle section z limits.
    """

    _default_title = 'Energy vs. z.'

    def __init__(
            self,
            piece,
            e_lims,
            tube_e_lims,
            title=None,
            z_lims=None,
            tube_z_lims=None,
            middle_z_lims=None
    ):
        super().__init__(piece, title)

        self.__e_lims = e_lims
        self.__tube_e_lims = tube_e_lims

        # Attributes with defaults.
        self.__z_lims = z_lims or z_lims[(8, 4)]
        self.__tube_z_lims = tube_z_lims or _default_tube_z_lims
        self.__middle_z_lims = middle_z_lims or _default_middle_z_lims

        self.__bins = _make_bins(
            self.__z_lims[0], self.__z_lims[1], self.bin_density
        )
        self.__bin_mids = _make_bin_midpoints(self.__bins)

    def plot_single(self, i=None, save=True, energy_label=None):
        """
        Plot a single event.

        :param energy_label: Energy deposit in middle tube sections,
            displayed in a text label on the plot.
        """
        fig, ax, ax_middle = self._make_fig_and_axes()
        fig.suptitle(self.title)
        self.__label_middle(ax_middle, energy_label)

        tube_kwargs = {
            'linewidth': _linewidth, 'color': 'purple', 'label': 'Tube Cals'
        }
        plate_kwargs = {
            'linewidth': _linewidth, 'color': 'blue', 'label': 'Plate Cals'
        }
        self.__split_plot(
            x=self.__bin_mids,
            ys=(self._to_density(self.get(i)),),
            plot_fn=ax_middle.plot,
            main_ax=ax,
            kwargs1=tube_kwargs,
            kwargs2=plate_kwargs
        )

        self.__add_legend(fig, ax, ax_middle)
        if save:
            self.save_fig(fig, file_suffix=None)

    def plot_means(self, energy_label=None):
        fig, ax, ax_middle = self._make_fig_and_axes()
        fig.suptitle(self.title + ' Averaged.')
        self.__label_middle(ax_middle, energy_label)

        means = self._to_density(numpy.mean(self.resultss, axis=0))
        stds = self._to_density(numpy.std(self.resultss, axis=0))

        tube_means_kwargs = {
            'linewidth': _linewidth,
            'color': 'purple',
            'label': 'Tube Cals Average'
        }
        plate_means_kwargs = {
            'linewidth': _linewidth,
            'color': 'blue',
            'label': 'Plate Cals Average'
        }
        tube_stds_kwargs = tube_means_kwargs.copy()
        plate_stds_kwargs = plate_means_kwargs.copy()
        tube_stds_kwargs.update({'alpha': 0.3, 'label': 'Standard Deviation'})
        plate_stds_kwargs.update({'alpha': 0.3, 'label': 'Standard Deviation'})

        self.__split_plot(  # Plot means.
            x=self.__bin_mids,
            ys=(means,),
            plot_fn=ax_middle.plot,
            main_ax=ax,
            kwargs1=tube_means_kwargs,
            kwargs2=plate_means_kwargs
        )
        self.__split_plot(  # Plot standard deviations.
            x=self.__bin_mids,
            ys=(means + stds, means - stds),
            plot_fn=ax_middle.fill_between,
            main_ax=ax,
            kwargs1=tube_stds_kwargs,
            kwargs2=plate_stds_kwargs
        )

        self.__add_legend(fig, ax, ax_middle)
        self.save_fig(fig, file_suffix=None)

    def plot_multi(self, energy_label=None):
        """Plot all of the events on one graph at once."""
        fig, ax, ax_middle = self._make_fig_and_axes()
        num_events = len(self.resultss)
        fig.suptitle(self.title + f' {num_events} events.')
        self.__label_middle(ax_middle, energy_label)

        kwargs = {'linewidth': _linewidth, 'alpha': 0.2}
        colors = itertools.cycle((
            'red', 'orange', 'yellow', 'pink', 'purple', 'brown', 'black'
        ))
        for _, sums in self.resultss.iterrows():
            kwargs['color'] = next(colors)
            self.__split_plot(
                x=self.__bin_mids,
                ys=(self._to_density(sums),),
                plot_fn=ax_middle.plot,
                main_ax=ax,
                kwargs1=kwargs,
                kwargs2=kwargs
            )

        self.save_fig(fig, file_suffix=f'{num_events}events')

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

    def _make_fig_and_axes(self):
        fig, ax = pyplot.subplots()
        ax_middle = ax.twinx()

        ax.set_xlabel(f'z ({self.length_units})')
        ax.set_ylabel(
            'E Dep Density - Plate Cals'
            f' ({self.energy_units} / {self.length_units})'
        )
        ax_middle.set_ylabel(
            'E Dep Density - Tube Cals'
            f' ({self.energy_units} / {self.length_units})'
        )

        # Indicate the middle tube sections with vertical bars and axis
        # ticks.
        for lim in self.__middle_z_lims:
            ax_middle.axvline(
                lim, linestyle='--', linewidth=1, color='gray', zorder=0
            )
        ax_middle.set_xticks(self.__middle_z_lims, minor=True)
        ax_middle.set_xticklabels(self.__middle_z_lims, minor=True)

        # for axes in (ax, ax_middle):
        #     axes.spines['right'].set_visible(False)
        #     axes.spines['top'].set_visible(False)

        for axes, e_limits in zip(
                (ax, ax_middle), (self.__e_lims, self.__tube_e_lims)
        ):
            axes.set_xlim(self.__z_lims)
            axes.set_ylim(e_limits)

        return fig, ax, ax_middle

    def __split_plot(
            self, x, ys, plot_fn, main_ax, kwargs1=None, kwargs2=None
    ):
        """
        Split up data and plot it on two axes.

        The data with x values inside split_lims (inclusive) is plotted
        using `plot_fn` normally, and the rest is scaled to `main_ax`.

        :param x: x values.
        :type x: 1d array
        :param ys: All the y values.
        :type ys: 2d array
        :param plot_fn: Function that plots data to the middle axes.
            The middle axes are displayed on top, so all plotting uses
            this function to avoid overlapping axes.
        :type plot_fn: function
        :param main_ax: The main axes, which lie underneath the middle
            axes.
        :type ax2_plot_func: matplotlib.pyplot.axes
        :param kwargs1: Keyword arguments passed to the middle plot.
        :type kwargs1: dict or None
        :param kwargs2: Keyword arguments passed to the full plot.
        :type kwargs2: dict or None
        :return: The two plot results.
        :rtype: (handle, handle)
        """
        kwargs1 = kwargs1 or None
        kwargs2 = kwargs2 or None

        inside_indices = _range2indices(x, self.__tube_z_lims)
        outside_indices = numpy.logical_not(inside_indices)

        plots = plot_fn(
            x[outside_indices],
            *(y[outside_indices] for y in ys),
            transform=main_ax.transData,
            **kwargs2
        )
        plots_middle = plot_fn(
            x[inside_indices],
            *(y[inside_indices] for y in ys),
            **kwargs1
        )
        return plots, plots_middle

    def __label_middle(self, ax, energy_label):
        """Label the energy deposit in the middle tube sections."""
        if energy_label:
            ax.annotate(
                energy_label,
                bbox={
                    'boxstyle': 'round4',
                    'facecolor':
                        pyplot.style.library['seaborn']['axes.facecolor'],
                    'edgecolor': 'silver',
                    'alpha': 0.7
                },
                xy=(0, 0),
                xycoords='data',
                xytext=(0, 25),
                textcoords='offset points',
                horizontalalignment='center',
                verticalalignment='center'
            )

    @staticmethod
    def __add_legend(fig, ax, ax_middle):
        """Place a legend at the upper-right corner of the axes"""
        # Combine legend handles and labels from both axes.
        handles, labels = (val1 + val2 for val1, val2 in zip(
            ax.get_legend_handles_labels(),
            ax_middle.get_legend_handles_labels()
        ))

        legend = ax_middle.legend(
            handles,
            labels,
            loc='upper right',
            ncol=2,
            fontsize='small',

            # Place legend at the corner of the axes rather.
            # `figure.legend` places the legend at the corner of the
            # figure instead.
            # bbox_to_anchor=(1, 1),
            # bbox_transform=ax.transAxes
        )
        legend.set_zorder(1)


class EnergyVsXY(Histogram):
    """Histogram of energy deposit vs. x and y."""

    _default_title = 'Energy vs. x and y.'

    def __init__(self, piece, title=None, xy_lims=None, z_lims=None):
        super().__init__(piece, title)
        self.resultss = []  # Ust a list instead of a DataFrame.

        # Attributes with defaults.
        self.__xy_lims = xy_lims or _default_xy_lims
        self.__z_lims = z_lims or _default_xy_hist_z_lims

        self.__binss = tuple(
            _make_bins(limits[0], limits[1], self.bin_density)
            for limits in self.__xy_lims
        )
        self.__binss_mesh = numpy.meshgrid(*self.__binss)
        # self.__bin_midss = tuple(
        #     _make_bin_midpoints(bins) for bins in self.__binss
        # )

    def plot_single(self, i=None, save=True):
        fig, ax = self._make_fig_and_axes()
        fig.suptitle(self.title)

        ax.pcolormesh(*self.__binss_mesh, self._to_density(self.get(i)))

        if save:
            self.save_fig(fig, file_suffix='xy')

    def plot_means(self):
        fig, ax = self._make_fig_and_axes()
        fig.suptitle(self.title + ' Averaged.')

        ax.pcolormesh(
            *self.__binss_mesh,
            self._to_density(numpy.mean(self.resultss, axis=0))
        )

        self.save_fig(fig, file_suffix='xy')

    def get(self, i=None):
        """
        Overrides the base to get from a list instead of a DataFrame.
        """
        if i is None:
            self._assert_single_entry()
            return self.resultss[0]

        return self.resultss[i]

    def add_results(self, sums):
        """
        Overrides the base to append to a list instead of a DataFrame.
        """
        self._check_results(sums)
        self.resultss.append(sums)

    def _data2results(self, data):
        sliced = data[numpy.logical_and(
            self.__z_lims[0] <= data.z, data.z <= self.__z_lims[1]
        )]
        # Transpose the sums to match the xy/Cartesian indexing used by
        # `numpy.meshgrid` and Matplotlib's `pcolormesh`.
        return numpy.histogram2d(
            sliced.x, sliced.y, bins=self.__binss, weights=sliced.energy_deposit
        )[0].T

    def _check_results(self, sums):
        """Check that `sums` is the correct shape for the bins."""
        assert numpy.all(numpy.equal(
            sums.shape,
            tuple(len(bins) - 1 for bins in self.__binss)
        ))

    def _make_fig_and_axes(self):
        fig, ax = pyplot.subplots()

        ax.set_xlabel(f'x ({self.length_units})')
        ax.set_ylabel(f'y ({self.length_units})')

        # No need to set x limits since the axes will be equal.
        # ax.set_xlim(self.__xy_lims[0])
        ax.set_ylim(self.__xy_lims[1])

        ax.set_aspect('equal')

        return fig, ax


def save_dataframe(dataframe, path):
    """
    Write a DataFrame to a csv file with nicely spaced columns and a
    header, and also write it to a second file compatible with Office
    365 (since it seems to prohibit the use of csv's).

    :param dataframe: The DataFrame to save.
    :param path: Path to save the files at. The extension of the
        filename at `path` is ignored and replaced with 'csv' and
        'xlsx'.
    :return: The paths to the saved files, (csv_path, excel_path).
    """
    tail, ext = os.path.splitext(path)

    csv_path = tail + '.csv'
    header = (
        'FCal and SCal analysis output.'
        f' {datetime.datetime.now().ctime()}'
    )
    with open(csv_path, 'w') as file:
        file.write('\n'.join((
            header,
            dataframe.to_string(index=False),
            ''
        )))

    excel_path = tail + '.xlsx'
    dataframe.to_excel(excel_path, index=False)

    return csv_path, excel_path


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


def _ordered_append(dataframe, other):
    """
    Append `other` to `dataframe` while keeping the order of the columns
    in both.

    The `append` method that comes with pandas Dataframes alphabetizes
    the order of the columns, which may not always be desired. This
    function provides an alternative.

    :param dataframe: The DataFrame to append to.
    :type dataframe: pandas.DataFrame
    :param other: The data to append.
    :type other: pandas.DataFrame or pandas.Series/dict-like object,
        or list of these
    :return: New DataFrame with appended data.
    :rtype: pandas.DataFrame
    """
    # The new column order. New columns from `other` are inserted in
    # front of the original `dataframe` columns.
    columns = [
                  col for col in other.keys() if col not in dataframe.keys()
              ] + list(dataframe.keys())

    # pandas' append alphabetizes the columns here. Indexing the result
    # with `columns` sets our own column order.
    return dataframe.append(other, ignore_index=True)[columns]
