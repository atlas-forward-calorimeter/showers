"""Analyzing data in "pieces." Events, Runs, etc.

TODO: Consider switching to a 16x9 aspect ratio.
TODO: Remove titles completely?

Changelog:
    append_mean_and_uncertainties
    combine_results
    turned off histograms
"""

import abc
import math
import os
import re
import textwrap
import traceback
import warnings

import numpy
import pandas

from . import calc

# PEEK box z-limits.
_PEEK_z_lims = (-43.5 / 2, 43.5 / 2)

# Position of the center of the "top right" calorimeter.
_tube_separation = 7.5
_single_cal_position = (
    3 / 4 * _tube_separation,
    math.sqrt(3) / 4 * _tube_separation,
    0
)

# Filename to use when saving combined `Numbers` results from several
# runs.
_numbers_filename = 'all-runs'

# Skip this many lines at the start when reading hits data files.
_skiprows = 1

# Regex pattern for numbers represented in filenames. An 'm' at the
# beginning means a minus/negative number. A 'p' stands for a decimal
# point.
_num_pattern = r'm?[0-9]*p?[0-9]+'


class Piece(abc.ABC):
    """
    A "piece" of analysis.

    TODO: Document this more.

        Attributes:
            name: A nice name to go by.

            hists: Contains all of the Histograms.

            z_title_width: Wrap energy-z plot titles to fit this many
                characters.

            xy_title_width: Same as `z_title_width` for energy-xy plots.
    """

    z_title_width = 58
    xy_title_width = 33

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
            if self.out_dir is None:
                self.out_dir = self.parent.out_subdir

        if self.out_dir:
            # Create out_dir in the filesystem and make a subdirectory
            # name.
            os.makedirs(self.out_dir, exist_ok=True)
            self.out_subdir = os.path.join(self.out_dir, self.name)
        else:
            self.out_subdir = None

        self._update_info()

        self.hists = {}

        z_title, xy_title = self._get_titles()
        self.hists['energy_vs_z'] = calc.EnergyVsZ(self, z_title)
        self.hists['energy_vs_xy'] = calc.EnergyVsXY(self, xy_title)

        self.numbers = calc.Numbers(self)

    def __repr__(self):
        return f'<{self.__class__.__name__ }: {self.name} @ {self.input_path}>'

    def _get_titles(self):
        """Get histogram titles from the Piece info."""
        title = 'E dep density.{}'

        if self.info.get('offset_beam'):
            title += ' Offset $e^-.$'

        tubes_angle = self.info.get('tubes_angle')
        if tubes_angle:
            title += f' ${tubes_angle}\\degree$ rotation.'

        plates = self.info.get('plates')
        if plates:
            title += f' {plates[0]} by {plates[1]}.'

        tube_shift = self.info.get('tube_shift')
        if tube_shift:
            title += f' {tube_shift} mm tube shift.'

        beam_shift = self.info.get('beam_shift')
        if beam_shift:
            title += f' {beam_shift} mm beam shift.'

        incident_energy = self.info.get('incident_energy')
        if incident_energy == '350gev':
            title += ' 350 GeV $e^-.$'
        elif incident_energy == '200gev':
            title += ' 200 GeV $e^-.$'

        if self.info.get('no_cryo'):
            title += ' No cryo.'

        # Add energy density descriptions and wrap the titles to meet
        # the `title_width`.
        z_title, xy_title = (
            '\n'.join(textwrap.wrap(title.format(description), width=width))
            for description, width in (
                ('', self.z_title_width),
                (' Middle section.', self.xy_title_width)
            )
        )

        return z_title, xy_title

    def _update_info(self):
        """Update this Piece's info (with access to all the defaults
        that are set in __init__). Implemented in the concrete Pieces.
        """

    @staticmethod
    def _energy_label(e_dep, std=None):
        if std:
            return f"{int(round(e_dep))}$\\pm${int(round(std))} MeV"
        return f"{int(round(e_dep))} MeV"

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


class Event(Piece):
    """
    A single event.

        Attributes:
            treated_hits: Hits data with the data in the PEEK box
             rotated so the tube cals are along the z-axis again.

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
        :type out_dir: str or False or None
        :param out_text: Path to output text. Defaults to
            out_dir/analysis.txt.
        :type out_text: str or None
        :param info: Specific parameters for this event.
        :type info: dict or None
        """
        super().__init__(hits_path, parent, out_dir, info)

        self.treated_hits = None

        # For tagging Numbers entries.
        self.__tags = {'event': self.info['event']}
        if self.parent:
            self.__tags['run'] = self.parent.name

        self.go()

    def go(self):
        """Read in data, make plots and do calculations."""
        hits = pandas.read_csv(self.input_path, skiprows=_skiprows)
        self.treated_hits = self.__treat_hits(hits)

        # for histogram in self.hists.values():
        #     histogram.add_data(self.treated_hits)
        self.numbers.add_data(self.treated_hits, self.__tags)

        if self.info.get('offset_beam'):
            e_dep = self.numbers.get()['top_right']
        else:
            e_dep = self.numbers.get()['middle_e_dep']

        # self.hists['energy_vs_z'].plot_single(
        #     energy_label=self._energy_label(e_dep)
        # )
        # self.hists['energy_vs_xy'].plot_single()

    def _update_info(self):
        event_number = int(self.name.split('-')[-1])
        self.info['event'] = event_number

    def _get_titles(self):
        z_title, xy_title = (
            f"{title} Event {self.info['event']}."
            for title in super()._get_titles()
        )
        return z_title, xy_title

    def _check_input_path(self, input_path):
        assert os.path.isfile(input_path), \
            f"The event data path {input_path} isn't a file."

    def _input_path2name(self, input_path):
        return _file2name(input_path)

    def __treat_hits(self, hits):
        """Reorient the data to reverse rotations and shifts of the tube
         cals.

        Rotates the PEEK hit positions about the center of the singled
        out tube until the tubes are aligned with the z-axis (as they
        were before the rotation tests), then cuts off any rotated
        data that ends up outside of the PEEK box limits.

        Reverses shifts of the tube cal(s) (not the offset of the beam).
        """
        offset_beam = self.info.get('offset_beam')
        tubes_angle = self.info.get('tubes_angle')
        tube_shift = self.info.get('tube_shift')
        if not offset_beam and not tubes_angle and not tube_shift:
            return hits

        in_PEEK = hits.z.between(*_PEEK_z_lims)

        PEEK_hits = hits[in_PEEK].copy()
        plate_hits = hits[~in_PEEK]
        PEEK_positions = PEEK_hits[['x', 'y', 'z']]

        if tube_shift:
            # Reverse the tube cal(s) shift.
            PEEK_positions -= (tube_shift, 0, 0)

        if tubes_angle:
            if offset_beam:
                # Shift so the center of the singled out tube cal is at the
                # origin.
                PEEK_positions -= _single_cal_position

            # Reverse the rotation about the y-axis.
            angle_radians = math.pi / 180 * tubes_angle
            reverseYrotation = numpy.array([
                [math.cos(-angle_radians),    0, math.sin(-angle_radians)],
                [0,                         1, 0                     ],
                [-math.sin(-angle_radians),   0, math.cos(-angle_radians)]
            ])
            PEEK_positions = (reverseYrotation @ PEEK_positions.T).T

            if offset_beam:
                # Shift the center of the singled out tube cal back to its
                # original position.
                PEEK_positions += _single_cal_position

        PEEK_hits[['x', 'y', 'z']] = PEEK_positions

        if tubes_angle:
            # Cut off the `PEEK_hits` that were rotated out of the PEEK
            # box range.
            PEEK_hits = PEEK_hits[PEEK_hits.z.between(*_PEEK_z_lims)]

        if offset_beam:
            # Cut out hits that aren't in the top right tube cal.
            # TODO: Organize the 7.5 / 4 offset!
            offset = 7.5 / 4
            PEEK_hits = PEEK_hits[PEEK_hits.y > 0]
            PEEK_hits = PEEK_hits[PEEK_hits.x > offset]

        treated_hits = pandas.concat(
            (PEEK_hits, plate_hits), ignore_index=True
        )

        return treated_hits


class Run(Piece):
    """A collection of events with a given setup."""

    # For "backwards compatibility" with the old data folder names.
    __old_plate_names = {'8-4': (8, 4), '7-5': (7, 5), '6-6': (6, 6)}

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

        means, stds, sems = self.numbers.append_mean_and_uncertainties(
            mean_tags={'run': self.name, 'event': 'mean'},
            std_tags={'run': self.name, 'event': 'std_dev'},
            sem_tags={'run': self.name, 'event': 'std_err'}
        )
        self.numbers.save()

        e_dep_key = (
            'top_right' if self.info.get('offset_beam') else 'middle_e_dep'
        )
        energy_label = self._energy_label(
            means[e_dep_key], stds[e_dep_key]
        )
        # self.hists['energy_vs_z'].plot_means(energy_label=energy_label)
        # self.hists['energy_vs_z'].plot_multi(energy_label=energy_label)
        # self.hists['energy_vs_xy'].plot_means()

        # Print out that we're done.
        print(f'Run {self.name} is finished!')

    def _new_events(self):
        """
        Create new analyzed Events and update this Run's analysis
        results.
        """
        for hits_path in self._event_paths():
            event = Event(hits_path, parent=self, out_dir=False)

            # for hist_name, hist in event.hists.items():
            #     self.hists[hist_name].add_results(hist.get())

            self.numbers.add_results(event.numbers.get())

    def _event_paths(self):
        """Get the paths to the event data files."""
        root, dirs, files = next(os.walk(self.input_path))
        for filename in files:
            if 'hits' in filename:
                yield os.path.join(root, filename)

    def _update_info(self):
        # TODO: Test that this still works.
        self.name2info(self.name, self.info)

    def _check_input_path(self, input_path):
        assert os.path.isdir(input_path), (
            f"The run directory {input_path} isn't actually a"
            f" directory."
        )

    def _input_path2name(self, input_path):
        return _dir2name(input_path)

    @staticmethod
    def name2info(name, info=None):
        """
        Parse the folder `name` and create or update the `info`
        dictionary (without replacing existing entries in `info`).

        If any of 8-4', '7-5', and '6-6' are in the filename, the first
        one of these, in this order, is taken as the plate arrangement.
        """
        if info is None:
            info = {}
        lowercase_name = name.lower()

        if 'incident_energy' not in info:
            if '350gev' in lowercase_name:
                info['incident_energy'] = '350gev'
            else:
                info['incident_energy'] = '200gev'
        if 'plates' not in info:
            # TODO: Double check this (quickly).
            plates_strings = re.findall(r'[0-9]+by[0-9]+', lowercase_name)
            assert len(plates_strings) < 2, \
                "Matched two or more plate arrangements."

            if plates_strings:
                plate_string = plates_strings[0]
                info['plates'] = tuple(
                    int(plate) for plate in plate_string.split('by')
                )
            else:
                for string, plates in Run.__old_plate_names.items():
                    if string in lowercase_name:
                        info['plates'] = plates
                        break
        if 'no_cryo' not in info and 'nocryo' in lowercase_name:
            info['no_cryo'] = True
        if 'offset_beam' not in info and 'offset' in lowercase_name:
            info['offset_beam'] = True
        # Rotation angle of the tube cals about the y-axis.
        if 'tubes_angle' not in info:
            angle_strings = re.findall(_num_pattern + 'deg', lowercase_name)
            if angle_strings:
                info['tubes_angle'] = _parse_number(angle_strings[0])
        if 'tube_shift' not in info:
            tube_shift_strings = re.findall(
                'tube-?' + _num_pattern + 'mm', lowercase_name
            )
            if tube_shift_strings:
                info['tube_shift'] = _parse_number(tube_shift_strings[0])
        if 'beam_shift' not in info:
            beam_shift_strings = re.findall(
                r'beam-?' + _num_pattern + 'mm', lowercase_name
            )
            if beam_shift_strings:
                info['beam_shift'] = _parse_number(beam_shift_strings[0])

        return info


def go(out_dir, *run_dirs, info=None):
    """
    Analyze a bunch of runs and output their results in an organized
    folder.

    :param run_out_dir: Where to save all the results, if at all.
    :type run_out_dir: str
    :param run_dirs: Locations of the run data directories.
    :type run_dirs: [str]
    :return:
    :rtype:
    """
    runs = []
    for run_dir in run_dirs:
        try:
            run_out_dir = os.path.join(out_dir, _dir2name(run_dir))
            runs.append(Run(
                events_path=run_dir, out_dir=run_out_dir, info=info
            ))
        except Exception as e:
            msg = f"Couldn't analyze the directory {run_dir}."
            warnings.warn(msg)
            print("Traceback (most recent call last):")
            traceback.print_tb(e.__traceback__)
            print(e)

    assert runs, "Couldn't successfully analyze any runs!"

    # if out_dir:
    #     calc.save_dataframe(
    #         dataframe=pandas.concat(
    #             (run.numbers.resultss for run in runs),
    #             ignore_index=True
    #         ),
    #         path=os.path.join(out_dir, _numbers_filename)
    #     )

    return runs


def combine_results(resultsss, out_file=None):
    """
    Combine analysis `resultss` from multiple runs, keeping the means
    and uncertainties from each run.

    :param resultsss: Sequence of `resultss` (2D DataFrame of results).
    :param out_file: Where to save the compiled means and uncertainties.
    :return: DataFrame of compiled means and uncertainties.
    """
    # Combine pandas Series as rows.
    combined = pandas.concat(
        (_clean_resultss(resultss) for resultss in resultsss),
        axis=1
    ).T

    # Fill in defaults.
    combined.fillna({
            'plates': '(8, 4)',
            'tubes_angle': 0,
            'beam_shift': 0,
            'tube_shift': 0
        }, inplace=True)

    if out_file:
        calc.save_dataframe(combined, out_file)
    return combined


def _clean_resultss(resultss, info=None):
    """Keep the means and uncertainties of a run's `resultss`, collapse
    them to one row, and add entries for the run's `info`.

    :type resultss: pandas.DataFrame
    :type info: dict
    :rtype: pandas.Series
    """
    # The row `tags` for the means and uncertainties.
    stats_rows = ('mean', 'std_dev', 'std_err')

    # Keep the means and uncertainties of these results. The means and
    # uncertainties will be combined and put on a single row.
    numbers_columns = ['middle_e_dep', 'tube_e_dep', 'full_e_dep']

    run_name = resultss['run'].iloc[0]

    if not info:
        info = Run.name2info(run_name)

    # Keep only the mean, standard deviation, and standard error rows.
    keep = resultss['event'].apply(
        lambda event: event in stats_rows
    )
    kept = resultss[keep]

    # Collapse to a single row (pandas Series). Add prefixes.
    clean = pandas.Series({'run_name': run_name})
    clean = clean.append(pandas.Series(info))
    for stat in stats_rows:
        clean = clean.append(
            kept[kept['event'] == stat][numbers_columns].squeeze().add_prefix(
                stat + '_')
        )

    return clean


def _file2name(file_path):
    """
    Return the name of the file at `file_path` without its extension.

    :param file_path: File path.
    :type file_path: str
    :return: Nice name.
    :rtype: str
    """
    tail, head = os.path.split(file_path)
    assert head != '', "Is this a directory instead of a file_path?"

    return head.split('.')[0]


def _dir2name(dir_path):
    """
    Return the name of the directory at `dir_path`.

    :param dir_path: Directory path.
    :type dir_path: str
    :return: Nice name.
    :rtype: str
    """
    tail, head = os.path.split(dir_path)
    if head == '':
        tail, head = os.path.split(tail)

    return head


def _parse_number(string, clean=True):
    """Parse and return a number from a string representation. This
    convention is used in the data filenames.

    :param string: String to parse and turn into a number.
    :param clean: If `True`, pick out the part of `string` that matches
        the number pattern.
    :return: The number obtained from `string`.
    """
    if clean:
        matches = re.findall(_num_pattern, string)
        assert matches, \
            f"Couldn't clean up this unrecognized string: {string}"
        assert len(matches) == 1, \
            f"Can't parse the multiple number patterns in {string}"
        string = matches[0]
    assert re.fullmatch(_num_pattern, string), \
        f"Unrecognized number format: {string}"

    # Interpret the leading 'm' as a minus sign and remove it.
    if string.startswith('m'):
        sign = -1
        string = string[1:]
    else:
        sign = 1

    # Replace 'p' with decimal points and convert the string to a float.
    number = sign * float(string.replace('p', '.'))

    if number.is_integer():
        number = int(number)
    return number


def _weak_update(dict1, dict2):
    """
    Update dict1 with dict2 but with reversed priority. The existing
    keys and values in dict1 are kept.
    """
    new_dict1 = dict2.copy()
    new_dict1.update(dict1)
    dict1.update(new_dict1)
