"""Analyze data from the command line.

This module is run when the package itself is executed in the command
line, e.g., with the command:
    $ python ./showers/
"""
from analysis import pieces

import argparse
from matplotlib import pyplot


def _showers():
    info = {}
    if args.incident_energy:
        info['incident_energy'] = _parse_incident_energy(args.incident_energy)

    runs = pieces.go(args.run_dirs, args.out, info=info, run_klass=pieces.Run)
    for run in runs:
        print(f'Run {run.name}:')
        print(run.numbers.resultss.tail())


def _foil():
    from analysis import foil
    multirun = foil.MultiRun(args.run_dirs, out_dir=args.out)
    pyplot.show()


def _parse_incident_energy(incident_energy):
    """Turn command line arguments into proper `info` items for
    the `incident_energy`.
    """
    if not incident_energy:
        return None

    incident_energy = incident_energy.lower()
    if incident_energy in ('350', '350gev'):
        return '350gev'
    elif incident_energy in ('200', '200gev'):
        return '200gev'
    else:
        return None


_parser = argparse.ArgumentParser()

styles = _parser.add_mutually_exclusive_group()
styles.add_argument(
    '-s',
    '--showers', action='store_const', const=_showers, dest='go',
    default=_showers,
    help="""Default mode. Perform the 'showers' analysis, which totals
    and plots energy deposits from high energy incident betas.
    """
)
styles.add_argument(
    '-f', '--foil', action='store_const', const=_foil, dest='go',
    help="""Perform the 'foil' analysis, which plots energy deposits 
    from foil sources with varying radii.
    """
)

_parser.add_argument(
    'run_dirs',
    nargs='*',
    help="""Directories containing hits data from single runs. At least
    one directory is required.
    """
)
_parser.add_argument(
    '-o', '--out',
    help="""Output results to this directory. The directory will be 
    created if it doesn't exist. If it does exist, files inside 
    may be overwritten.
    """
)
_parser.add_argument(
    '-e', '--incident_energy',
    help="""Sets the energy of the incident electrons. '200', '200GeV', 
    '350', and '350GeV' are currently recognized. Case is ignored.
    """
)

args = _parser.parse_args()

if args.run_dirs:
    args.go()
else:
    _parser.print_usage()
    print('At least one run directory is required. Add -h or --help to '
          'display help.')
