"""The `Piece` and `FCalPiece` base classes."""

import os
import warnings

from .helpers import makeBins, binMidpoints

class Piece:
    """A "piece" of analysis.

    An analysis procedure can be thought of as being made up of pieces.
    Each piece can be made up of smaller pieces and make up bigger ones.
    For example, a collection of runs might be one big piece, which is
    made up of runs that are smaller pieces, which are each made up of
    events, the smallest piece.

    `smallerPieces(self)`:

    If this piece is made up of smaller ones, this function should be
    implemented in derived classes to give them. It should take no
    parameters other than `self`.

    return: An iterator over `Piece` objects of the next lowest 
    level of analysis.

    `self` should be assigned as the `parent` of the smaller piece
    in this method, if the smaller piece needs it.

    `analyze(self)`:

    Implement in derived classes to carry out analysis on the level of
    this piece.
        
    If this piece is made of smaller ones, `go` will have the smaller 
    pieces do their analysis first, then this method will finally be 
    called.

    """
    smallerPieces = None
    analyze = None

    def __init__(self, name=None, parent=None):
        """Initialize everything before looping through smaller pieces.

          name: Just to have a nicer sounding name to go by.
        parent: The next larger piece, if there is one. If not, this
                piece will be a "top level" piece.
        
        Note: Concrete derived classes will usually want to take some
        input and output parameters in their constructor. For example,
        an input might be a path or URL leading to some data, and an
        output could be similar.
        """
        self.name = name
        self.parent = parent

    def go(self):
        """Recursively analyze! 
        
        Analyzes each of the next smallest pieces before analyzing
        this piece itself. Call `go` on the biggest piece to analyze 
        everything in it.
        
        """
        if self.smallerPieces:
            # Loop over next smallest level.
            for piece in self.smallerPieces():
                piece.go()
        
        if self.analyze:
            # Analyze this level.
            self.analyze()
    
    def __repr__(self):
        # TODO: Finish this.
        if self.name:
            return f'<{type(self).__name__}: {self.name}>'
        else:
            return super().__repr__()


class FCalPiece(Piece):
    """A piece class with goodies specific to the FCal analysis.
    
    This is an abstract class, building extra functionality on top of 
    `Piece`.
    
    """
    # Geometry parameters.
    tubeZ = 35 / 2
    tubeMiddleZ = 8 / 2
    # Analysis parameters.
    fullZlims = (-50, 80)
    xyLims = 7
    binDensity = 10
    xyBinDensity = 10
    plotFileFormat = 'svg'

    # Analysis numbers derived from above numbers.
    fullBins = makeBins(*fullZlims, binDensity)
    fullBinMids = binMidpoints(fullBins)
    xyBins = makeBins(-xyLims, xyLims, binDensity)
    xyBinMids = binMidpoints(xyBins)

    def __init__(self, 
                 inputPath, 
                 outDirectory=None,
                 maxEvents=None,
                 parent=None):
        """inputPath: Path of input data. Could be a folder or a file.
        outDirectory: Write analysis output in here. Defaults to
                      `inputPath`/analysis for top level analysis. Set
                      to `False`, 0, etc. to write nothing to disk.
           maxEvents: The maximum number of events per run to analyze.
        
        """
        name = os.path.split(inputPath)[-1]
        super().__init__(name=name, parent=parent)

        if outDirectory is None and parent:
            # Propagate output paths from parent.
            outDirectory = parent.outDirectory
            outTextPath = parent.outTextPath
        else:
            if outDirectory is None:
                # By default, output to the directory
                # `inputPath`/analysis.
                outDirectory = os.path.join(inputPath, 'analysis')
            if outDirectory:
                # Write output text to `outDirectory`/analysis.txt.
                outTextPath = os.path.join(outDirectory, 'analysis.txt')
            else:
                outTextPath = None
        
        if maxEvents is None:
            if parent:
                # Propagate max number of events from parent.
                maxEvents = parent.maxEvents
        
        if not parent and not self.smallerPieces:
            warnings.warn(
                "Can't output results to file when analyzing an event by"
                " itself."
            )

        self.inputPath = inputPath
        self.outDirectory = outDirectory
        self.maxEvents = maxEvents
        self.outTextPath = outTextPath
