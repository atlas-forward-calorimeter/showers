'''
Analysis of the hits output by Geant4, using Python.

Written by Anson Kost with the help of Professor John Rutherfoord, May 2019.

'''

import sys
import pandas as pd
from matplotlib import pyplot as plt

df = pd.read_csv(sys.argv[1], names=[
    'energyDeposit', 'x', 'y', 'z', 
    'energy', 'momentumX', 'momentumY', 'momentumZ', 'trackID'
])

plt.figure()
df.z.hist(bins=500)

plt.figure()
df.z[df.z.abs() < 17.5].hist(bins=500)

plt.show()
