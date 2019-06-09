"""The `MultiRun` class."""

import os
import datetime
from matplotlib import pyplot as plt

from .fcalpiece import FCalPiece
from .run import Run
from .helpers import doublePrint, outTextPath

class MultiRun(FCalPiece):
    """A collection of runs.
    
    One level higher than a run. The top level.
    
    """
    def __init__(self, 
                 dataDirectory, 
                 outDirectory=None,
                 params=None,
                 name=None, 
                 maxEvents=None):
        """dataDirectory: The top level directory of data to analyze. 
                      Every subdirectory inside that contains only files
                      is treated as a directory containing data from one 
                      run.
        outDirectory: Write all analysis output in here. Defaults to
                      `dataDirectory`/analysis. Set to `False`, 0, etc.
                      (but not `None`) to write nothing to files.
            `params`: Analysis parameters. Geometry dimensions, 
                      constants, etc.
              `name`: Just to have a nicer sounding name to go by.             
           maxEvents: The maximum number of events per run to analyze.

        """
        super().__init__(dataDirectory, outDirectory)

        # Default parameters.
        self.params = {
            'full z limits': (-50, 80),
            'tube z': 35 / 2,
            'tube middle z': 8 / 2,
            'num events': 28,
            'plate xy': 40 / 2,
            'bin density': 10
        }
        if params:
            self.params.update(params)

        # Output to `dataDirectory`/analysis by default.
        if self.outDirectory is None:
            self.outDirectory = os.path.join(dataDirectory, 'analysis')

        self.dataDirectory = dataDirectory  # synonym for `inputPath`
        self.maxEvents = maxEvents

        self.initialize()
    
    def initialize(self):
        """Like a constructor, but for analysis and output."""
        if self.outDirectory:
            # Create output directory
            # (and overwrite any existing analysis).
            os.makedirs(self.outDirectory, exist_ok=True)
        
        # Analysis Header
        analysisHeader = ('FCal analysis output. '
                          + datetime.datetime.now().ctime())
        self.doublePrint(analysisHeader)
    
    def smallerPieces(self):
        """Get every subdirectory of `directory` that contains only
        files.
        
        These directories should contain run data. The output directory
        itself is skipped.
        
        """
        for root, dirs, files in os.walk(self.dataDirectory):
            # Skip output directory.
            if self.outDirectory:
                if os.path.samefile(root, self.outDirectory):
                    continue
            if not dirs:
                # `root` only contains files.
                yield Run(dataDirectory=root, 
                          outDirectory=self.outDirectory, 
                          params=self.params,
                          maxEvents=self.maxEvents,
                          parent=self)
    
    def analyze(self):
        # At the end of everything.
        print('Done. Safe to close.')
        plt.show()
