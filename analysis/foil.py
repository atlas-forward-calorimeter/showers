"""Patch for analyzing foil radius effects. See pieces.py and calc.py
for original documentation."""

from . import pieces, calc

import itertools
import time

import pandas
from matplotlib import pyplot

mean_sr90_energy = 568.514688585 / 1000  # average MeV per beta


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
    max_events = 100000

    def __init__(self, *args, show_plots=True, **kwargs):
        self.show_plots = show_plots

        super().__init__(*args, go=False, **kwargs)
        self.hists = {'energy_vs_z': EnergyVsZ(
                self, self._get_titles()[0], N=self.info.get('num_betas')
        )}
        self.numbers = None

        self.go()

    def go(self):
        # Read in event(s).
        switch = False
        multiprocess = True
        if switch:
            hits = pandas.read_csv('analysis/tests/10000-betas/hits-0.csv', skiprows=pieces.skiprows)
            hist = hits.groupby('track_id').apply(self.hists['energy_vs_z']._data2results)
            self.hists['energy_vs_z'].hist = hist
            self.hists['energy_vs_z'].add_data(hits)
        if not switch:
            if multiprocess:
                from multiprocessing import Pool
                with Pool(processes=28) as pool:
                    for sums in pool.map(
                        self._path2sums,
                        itertools.islice(self._event_paths(), self.max_events)
                    ):
                        if sums is not None:
                            self.hists['energy_vs_z'].add_results(sums)
            else:
                for hits_path in itertools.islice(
                        self._event_paths(), self.max_events
                ):
                    self.read_hits(hits_path)

        if self.show_plots:
            self.hists['energy_vs_z'].plot_means()

    @staticmethod
    def name2info(name, info=None):
        info = pieces.Run.name2info(name, info)
        if 'incident_energy' in info:
            del info['incident_energy']
        return info

    def _path2sums(self, hits_path):
        try:
            hits = pandas.read_csv(hits_path, skiprows=pieces.skiprows)
        except pandas.errors.EmptyDataError:
            print(f'{hits_path} is empty.')
            return None
        print(f'Read {hits_path}.')
        return self.hists['energy_vs_z']._data2results(hits)


class MultiRun(pieces.Piece):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.hists = {'energy_vs_z': EnergyVsZ(
            self,
            title='Energy deposit vs. z, averaged per decay, for different foil gap sizes.'
        )}

        # Make a single input_path iterable
        if isinstance(self.input_path, str):
            self.input_path = [self.input_path]

        self.go()

    def go(self):
        start = time.time()
        self._new_runs()
        stop = time.time()
        elapsed_mins = (stop - start) / 60
        print(
            f'{self.name} took {elapsed_mins} minutes to analyze '
            f'individual runs.'
        )
        self.hists['energy_vs_z'].pickle_save()
        self.hists['energy_vs_z'].plot_multi()

    def _new_runs(self):
        for path in self.input_path:
            run = Run(events_path=path, out_dir=self.out_dir, show_plots=False)
            for name, hist in self.hists.items():
                hist.add_mean_of_results(run.hists[name], info=run.info)

    def _check_input_path(self, input_path):
        pass

    def _input_path2name(self, input_path):
        return 'MultiRun'


class EnergyVsZ(calc.EnergyVsZ):

    figsize = (12.8, 7.2)

    def __init__(self, *args, N=None, **kwargs):
        self.N = N
        self.infos = []
        e_lims = (0, 6)
        tube_e_lims = e_lims
        super().__init__(
            *args,
            e_lims=e_lims,
            tube_e_lims=tube_e_lims,
            z_lims=(-20 / 2, 20 / 2),
            energy_units='KeV',
            bin_density=5,
            **kwargs
        )

    def plot_means(self, energy_label=None):
        self._assert_nonempty()

        fig, ax = self._make_fig_and_axes()
        fig.suptitle(self.title + ' Averaged.')
        self.__label_middle(ax, energy_label)

        if not self.N:
            self.N = self.resultss.shape[0]
        means = self._to_density(self.mean_per_particle())

        # stds = self._to_density(numpy.std(self.resultss, axis=0))
        # sems = stds / math.sqrt(self.resultss.shape[0])
        # sems = self._to_density(
        #     numpy.mean(numpy.sqrt(self.resultss * mean_sr90_energy), axis=0) / (N ** (2 / 2))
        # )
        # sems = means / N ** (1 / 4)
        # for val in self.stds:
        #     print(val)
        # M = self.hist.shape[0]
        # sems = 1000 * self._to_density(
        #     math.sqrt(M / N) * numpy.std(self.hist, axis=0) / math.sqrt(N)
        # )

        means_kwargs = {
            'linewidth': 2 * calc.linewidth,
            'color': 'red',
            'label': 'Average E Dep Density'
        }
        stds_kwargs = means_kwargs.copy()
        stds_kwargs.update({'alpha': 0.3, 'label': 'Standard Error'})

        ax.plot(self.__bin_mids, means, **means_kwargs)
        # ax.fill_between(
        #     self.__bin_mids, means + sems, means - sems, **stds_kwargs
        # )

        # ax.legend(loc='upper right')
        self.save_fig(fig, file_suffix='z')

        return ax

    def plot_multi(self):
        """Relies on `self.infos` and `self.resultss` corresponding to
        each other.
        """
        fig, ax = self._make_fig_and_axes()
        fig.suptitle(self.title)

        for (_, row), info, color in zip(
                self.resultss.iterrows(), self.infos, self.colors
        ):
            foil_gap = info.get('foil_gap')
            label = f'{foil_gap} mm' if foil_gap else None

            ax.plot(
                self.__bin_mids,
                self._to_density(row),
                label=label,
                color=color,
                linewidth=2 * calc.linewidth
            )

        # Sort legend by its labels.
        handles, labels = ax.get_legend_handles_labels()
        pairs = sorted(zip(handles, labels), key=lambda pair: pair[1])
        ax.legend(*zip(*pairs), title='Foil gap size:')

        self.save_fig(fig, file_suffix=None)

    def add_mean_of_results(self, other, info=None):
        """Add mean per particle results of another histogram. Use this
        when using `plot_multi`.

        :param other: Other `EnergyVsZ` histogram.
        :type other: EnergyVsZ
        """
        self.add_results(other.mean_per_particle())
        self.infos.append(info or {})

    def mean_per_particle(self, N=None):
        """Total results divided by number of particles."""
        if not N:
            if not self.N:
                raise RuntimeError(
                    "`N` (number of particles) must be defined to take "
                    "means per particle."
                )
            N = self.N
        return self.resultss.sum() / N

    @property
    def colors(self):
        return itertools.cycle((
            'red', 'blue', 'orange', 'purple', 'gray', 'brown', 'green',
            'black', 'darkturquoise'
        ))

    def _make_fig_and_axes(self):
        fig, ax = pyplot.subplots(figsize=self.figsize)

        ax.set_xlabel(f'z ({self.length_units})')
        ax.set_ylabel(
            'Energy Deposit Density'
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

        # ax.spines['right'].set_visible(False)
        # ax.spines['top'].set_visible(False)

        ax.set_xlim(self.__z_lims)
        ax.set_ylim(self.__e_lims)

        return fig, ax

    def _to_density(self, sums):
        """In KeV/mm rather than MeV/mm."""
        return 1000 * super()._to_density(sums)
