'''
Analysis of the hits output by Geant4, using Python.

Written by Anson Kost with the help of Professor John Rutherfoord, May 2019.

'''

import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# Read in data.
df = pd.read_csv(
    sys.argv[1],    # file path
    names = ('e',   # energy deposit   
             'x',   # x position
             'y',   # y position
             'z',   # z position
             'id')  # track ID
)

tubeZ = 35 / 2          # Half length of tubes in mm.
tubeMiddleZ = 8 / 2     # Half length of middle tube section in mm.

def analyze(name, dataFrame, kind='all'):
    '''
    Analyze one set of data.
    `binsPerMil`: Bins per millimeter for histgrams.
    `kind`: Kind of data - 'all' or just the 'tubes'.
    '''
    ## Energy vs. z histogram.
    # Plot's x and y boundaries.
    histZLimit = 70 if kind == 'all' else tubeZ
    histELimit = 35 if kind == 'all' else 2.25
    
    binsPerMil = 10  # Bins per millimeter for histograms.
    
    histBins = np.linspace(  # The bins.
        -histZLimit, histZLimit, int(binsPerMil * 2 * histZLimit + 1))
    
    # Energy vs. z histogram.
    plt.figure()
    plt.title(name)
    plt.ylim(0, histELimit)
    dataFrame.z.hist(weights=dataFrame.e, bins=histBins)
    
    ## Energy deposit sums.
    print(name)
    print(f'Total energy deposit: \n{dataFrame.e.sum()} MeV')
    
    if kind == 'tubes':
        print('Energy deposit in middle tubes: \n'
              f'{dataFrame[dataFrame.z.abs() < tubeMiddleZ].e.sum()}'
              ' MeV')


analyze('All', df)
analyze('Tubes', df[df.z.abs() < tubeZ], kind='tubes')

plt.show()
