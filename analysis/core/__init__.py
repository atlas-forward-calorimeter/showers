"""Analysis of the hits output by Geant4, using Python.

Each output file from Geant4 represents an event. Each folder represents a run.
"Internal units" are millimeters (mm) and megaelectronvolts (MeV).

Written by Anson Kost with the help of Professor John Rutherfoord, May 2019.

"""

from .utilities import *
from .analyses import *
