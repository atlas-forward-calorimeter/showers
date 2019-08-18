"""Handle arguments from the command line.

This module is run when the package itself is executed in the command
line, e.g., with the command:
    $ python ./showers/
"""

import argparse

from analysis import pieces

_parser = argparse.ArgumentParser()
_parser.add_argument(
    '-o',
    '--out',
    help="""Output results to this directory. The directory will be 
    created if it doesn't exist. If it does exist, files inside 
    may be overwritten.
    """
)
_parser.add_argument(
    '-e',
    '--incident_energy',
    help="""Sets the energy of the incident electron. '200', '200GeV', 
    '350', and '350GeV' are currently recognized. Case is ignored.
    """
)
_parser.add_argument('run_dirs', nargs='*')
args = _parser.parse_args()


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


_run_dirs = args.run_dirs or ['.']

_info = {}
if args.incident_energy:
    _info['incident_energy'] = _parse_incident_energy(args.incident_energy)

runs = pieces.go(args.out, *_run_dirs, info=_info)
for run in runs:
    print(run.numbers.resultss.tail())
