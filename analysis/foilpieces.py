"""Patch for analyzing foil radius effects. See pieces.py and calc.py
for original documentation."""

from . import pieces, calc

import math
import itertools

import numpy
import pandas
from matplotlib import pyplot


class Event(pieces.Event):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, go=False, **kwargs)
        print(f'Beginning event {self.name}.')

        self.empty_hits = False
        self.hists = {'energy_vs_z': EnergyVsZ(self)}
        self.numbers = None

        self.go()

    def go(self):
        try:
            hits = pandas.read_csv(self.input_path, skiprows=pieces.skiprows)
        except pandas.errors.EmptyDataError:
            self.empty_hits = True
            return

        for histogram in self.hists.values():
            histogram.add_data(hits)
        if self.numbers:
            self.numbers.add_data(hits, self.__tags)


class Run(pieces.Run):

    event_klass = Event
    max_events = 10000

    def __init__(self, *args, **kwargs):
        super().__init__(*args, go=False, **kwargs)
        self.hists = {'energy_vs_z': EnergyVsZ(self, self._get_titles()[0])}
        self.numbers = None
        self.go()

    def go(self):
        # self._new_events()

        switch = False

        if switch:
            hits = pandas.read_csv('analysis/tests/data3/hits-0.csv', skiprows=pieces.skiprows)
            self.hists['energy_vs_z'].add_data(hits)

        if not switch:
            for hits_path in itertools.islice(self._event_paths(), self.max_events):
                try:
                    hits = pandas.read_csv(hits_path, skiprows=pieces.skiprows)
                except pandas.errors.EmptyDataError:
                    print(f'{hits_path} is empty.')
                    continue
                print(f'Read {hits_path}.')
                self.hists['energy_vs_z'].add_data(hits)

        self.hists['energy_vs_z'].plot_means()

    @staticmethod
    def name2info(name, info=None):
        info = pieces.Run.name2info(name, info)
        if 'incident_energy' in info:
            del info['incident_energy']
        return info


class EnergyVsZ(calc.EnergyVsZ):
    def __init__(self, *args, **kwargs):
        e_lims = (0, 10)
        tube_e_lims = e_lims
        super().__init__(
            *args,
            e_lims=e_lims,
            tube_e_lims=tube_e_lims,
            z_lims=(-20 / 2, 20 / 2),
            energy_units='KeV',
            **kwargs
        )

    def plot_means(self, energy_label=None):
        self._assert_nonempty()
        fig, ax = self._make_fig_and_axes()
        fig.suptitle(self.title + ' Averaged.')
        self.__label_middle(ax, energy_label)

        # In KeV/mm rather than MeV/mm.
        N = 1
        means = 1000 * self._to_density(numpy.mean(self.resultss, axis=0)) / N
        stds = 1000 * self._to_density(numpy.std(self.resultss, axis=0))
        sems = stds / math.sqrt(self.resultss.shape[0])
        # sems = 1000 * self._to_density(numpy.mean(numpy.sqrt(self.resultss), axis=0) / (N ** (2 / 2)))

        means_kwargs = {
            'linewidth': calc.linewidth,
            'color': 'red',
            'label': 'Average E Dep Density'
        }
        stds_kwargs = means_kwargs.copy()
        stds_kwargs.update({'alpha': 0.3, 'label': 'Standard Error'})

        ax.plot(self.__bin_mids, means, **means_kwargs)
        ax.fill_between(
            self.__bin_mids, means + sems, means - sems, **stds_kwargs
        )

        ax.legend(loc='upper right')
        self.save_fig(fig, file_suffix='z')

    def _make_fig_and_axes(self):
        fig, ax = pyplot.subplots()
        # ax_middle = ax.twinx()

        ax.set_xlabel(f'z ({self.length_units})')
        ax.set_ylabel(
            'E Dep Density'
            f' ({self.energy_units} / {self.length_units})'
        )

        # Indicate the middle tube sections with vertical bars and axis
        # ticks.
        for lim in self.__middle_z_lims:
            ax.axvline(
                lim, linestyle='--', linewidth=1, color='gray', zorder=0
            )
        ax.set_xticks(self.__middle_z_lims, minor=True)
        ax.set_xticklabels(self.__middle_z_lims, minor=True)

        # for axes in (ax, ax_middle):
        #     axes.spines['right'].set_visible(False)
        #     axes.spines['top'].set_visible(False)

        # for axes, e_limits in zip(
        #         (ax, ax_middle), (self.__e_lims, self.__tube_e_lims)
        # ):
        ax.set_xlim(self.__z_lims)
        # ax.set_ylim(self.__e_lims)

        return fig, ax
