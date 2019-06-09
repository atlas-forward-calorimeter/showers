"""The `FCalPiece` base class."""

import os
from .piece import Piece

class FCalPiece(Piece):
    """A piece with goodies specific to the FCal analysis.
    
    This is an abstract class built on top of `Piece`.
    
    """
    numEvents = 28,
    
    fullZlims = (-50, 80)
    tubeZ = 35 / 2,
    tubeMiddleZ = 8 / 2,
    plateXY = 40 / 2,
    binDensity = 10

    def __init__(self, 
                 inputPath, 
                 outDirectory, 
                 isAtom=False, 
                 parent=None):
        """
        #TODO: Document this!
        """
        super().__init__(name=inputPath, isAtom=isAtom, parent=parent)
        self.inputPath = inputPath
        """
        self.outDirectory = outDirectory
        if outDirectory:
            self.outTextPath = os.path.join(outDirectory, 'analysis.txt')
        else:
            self.outTextPath = None
        
        # why does this need to be in every child?
        """

    def doublePrint(self, *objs, **kwargs):
        """Print a message and append it to a file at the same time!
        
        `objs`: Object(s) to print and write.
        `filePath`: File to write to. Doesn't write to a file by 
        default.
        
        `end` and keyword arguments are the same as in the built-in
        `print`.
        
        """
        print(*objs, **kwargs)
        if self.outTextPath:
            message = ' '.join(objs) + '\n'
            with open(self.outTextPath, 'a') as file:
                file.write(message)
