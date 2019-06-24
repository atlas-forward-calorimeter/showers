"""The `MultiRun` class."""

import os
import datetime
from matplotlib import pyplot as plt

from .piece import FCalPiece
from .run import Run
from .helpers import printAndWrite

class MultiRun(FCalPiece):
    """A collection of runs.
    
    One level higher than a run. The top level.
    
    """
    def __init__(self, 
                 dataDirectory, 
                 outDirectory=None,
                 maxEvents=None):
        """dataDirectory: The top level directory of data to analyze. 
                      Every subdirectory inside that contains only files
                      is treated as a directory containing data from one 
                      run.
           maxEvents: The maximum number of events per run to analyze.

        """
        super().__init__(
            dataDirectory, outDirectory=outDirectory, maxEvents=maxEvents)

        self.dataDirectory = dataDirectory  # same as `self.inputPath`
        self.maxEvents = maxEvents

        self.start()
    
    def start(self):
        """Like a constructor, but for analysis and output."""
        if self.outDirectory:
            # Create output directory
            # (and overwrite any existing analysis).
            os.makedirs(self.outDirectory, exist_ok=True)
            print('Creating out directory.')
        
        # Analysis Header
        analysisHeader = (
            f'FCal analysis output. {self.name} directory.'
            f' {datetime.datetime.now().ctime()}'
        )
        printAndWrite(analysisHeader, file=self.outTextPath)
    
    def smallerPieces(self):
        """Get every subdirectory of `directory` that contains only
        files.
        
        These directories should contain run data.
        
        """
        skipDirs = ['.git']
        if self.outDirectory:
            skipDirs.append(os.path.split(self.outDirectory)[-1])
        for root, dirs, files in os.walk(self.dataDirectory):
            # Skip directories in `skipDirs` 
            # (folders with the name ".git" or the same name as the 
            # output directory).
            dirs[:] = [d for d in dirs if d not in skipDirs]
            if not dirs:
                # `root` only contains files.
                yield Run(dataDirectory=root, parent=self)
    
    def analyze(self):
        # At the end of everything.
        printAndWrite('Done. Safe to close.', file=self.outTextPath)
        plt.show()
