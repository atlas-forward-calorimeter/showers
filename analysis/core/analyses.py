"""Different analysis functions with specific and general uses.

Written by Anson Kost with the help of Professor John Rutherfoord, May 2019.

"""

import types
import datetime

from .utilities import *

## Settings ##
fullZlimits = (-50, 80)
tubeZ = 35 / 2
tubeMiddleZ = 8 / 2
## ##

def analysis1(dataPath, 
              outPath=None, 
              outFileName='analysis.txt', 
              maxFiles=None):
    """Make energy-z histograms for events and for runs.

    `dataPath`:     Path to top level data folder containing secondary level 
                    folders of events.
    `outPath`:      Top folder to save analysis results in. 
                    Default is `dataPath`/analysis.
    `outFileName`:  Save calculation results to this filename inside `outPath`.
    'maxFiles':     Maximum # of events per run to analyze.
    
    This is the first analysis function made. It's simple and does all the
    looping on its own (with the help of the utility functions).

    """
    ## Histogram limits for full data. ##
    fullHistParams = {
        'zLimits': fullZlimits,    
        'binDensity': None,     # default 10 per mm
        'energyLimit': 38       # vertical scale
    }

    ## Histogram limits for tubes only. ##
    tubesHistParams = {
        'zLimits': (-tubeZ, tubeZ),    
        'binDensity': None,     # default 10 per mm
        'energyLimit': 3        # vertical scale
    }

    # Bins and midpoints based on above limits.
    fullBins = makeBins(*fullHistParams['zLimits'], 
                        fullHistParams['binDensity'])
    fullBinMidpoints = binMidpoints(fullBins)

    # Default output path.
    if outPath is None:
        outPath = os.path.join(dataPath, 'analysis')
    
    os.makedirs(outPath, exist_ok=True)  # Make output folders if needed.

    # Analysis output file.
    outFilePath = os.path.join(outPath, outFileName)

    # Analysis File Header
    analysisHeader = 'FCal python analysis output. ' \
                     + datetime.datetime.now().ctime() \
                     + '\n'
    doublePrint(analysisHeader, outFilePath)

    for folderName, folder in readFolder(dataPath, maxFiles):
        ## Run Level ##

        # Run level data and output paths.
        runDataPath = os.path.join(dataPath, folderName)
        runOutPath = os.path.join(outPath, folderName)

        # Skip the analysis output folder itself.
        if outPath == runDataPath:
            continue

        # Assert subfolders exist.
        assert(isinstance(folder, types.GeneratorType))

        os.makedirs(runOutPath, exist_ok=True)  # Make output folders if needed.

        # Run Header
        runHeader = f'\nBegin of run {folderName}.\n\n'
        doublePrint(runHeader, outFilePath)

        # "Initalize" run calculations.
        runFullEdep, runMiddleEdep = 0, 0

        # "Initialize" run histogram sums.
        runFullSums = np.zeros(len(fullBinMidpoints))

        # Multiple event histogram figure and formatting.
        plt.figure()
        plt.title(
            f'Histogram - Energy Deposit vs. z - Events - Run {folderName}')
        plt.xlabel('z')
        plt.ylabel('Energy Deposit Per Bin')

        for fileName, event in folder:
            ## Event Level ##

            # Assert there are files in 2nd level folders.
            assert(isinstance(event, pd.DataFrame))

            # Event level data and output paths.
            eventDataPath = os.path.join(runDataPath, fileName)
            eventOutPath = os.path.join(runOutPath, fileName)

            # Event Calculations
            fullEdep, middleEdep = calculate(event)
            # Cumulate run calculations.
            runFullEdep += fullEdep         
            runMiddleEdep += middleEdep

            # Write event calculations to analysis file.
            output = (
                f'Event {fileName}.\n'
                'Total Energy Deposit:\n'
                f'{fullEdep} MeV\n'
                'Middle Tube Energy Deposit:\n'
                f'{middleEdep} MeV\n'
                '-\n'
            )
            doublePrint(output, outFilePath)

            # Event Histogram(s)
            fullSums, _ = makeHist(event, fullBins)
            plt.plot(fullBinMidpoints, fullSums, lw=0.5)
            runFullSums += fullSums  # Cumulate run histogram.

        ## Run Level ##

        # Finish run histograms.
        plt.savefig(os.path.join(runOutPath, 'EventHist.svg'), format='svg')
        plt.figure()
        plt.title(f'Histogram - Energy Deposit vs. z - Sum - Run {folderName}')
        plt.xlabel('z')
        plt.ylabel('Energy Deposit Per Bin')
        plt.plot(fullBinMidpoints, runFullSums)
        plt.savefig(os.path.join(runOutPath, 'RunHist.svg'), format='svg')

        # Run calculations output.
        runOutput = (
            f'Run {folderName}.\n'
            'Run: Total Energy Deposit:\n'
             f'{runFullEdep} MeV\n'
            'Run: Middle Tube Energy Deposit:\n'
            f'{runMiddleEdep} MeV\n'
            '-\n'
        )

        # Run Footer
        runFooter = f'\nEnd of run {folderName}.\n\n'
        
        doublePrint(runOutput + runFooter, outFilePath)
    
    ## Analysis End ##
    print('Done. Safe to stop program.')
    plt.show()


def calculate(df):
    """Carry out calculations of important stuff."""
    # Energy deposit over all data.
    fullEnergyDeposit = df.energy_deposit.sum()
    # Energy deposit over middle tube electrode.
    middleEnergyDeposit = df.energy_deposit[df.z.abs() < tubeMiddleZ].sum()
    return fullEnergyDeposit, middleEnergyDeposit


def doublePrint(message, filePath):
    """Print to the terminal and append to a file at the same time!"""
    print(message)
    with open(filePath, 'a') as file:
        file.write(message)
