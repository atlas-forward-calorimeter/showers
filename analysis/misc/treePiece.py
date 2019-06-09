"""The `TreePiece` base class."""

class TreePiece(Piece):
    """Another base class that is just for analyses of data organized in
    directory trees.

    """
    def __init__(self, putIn, putOut, isAtom=False):
        """ putIn: The directory containing data to analyze. 
           putOut: Write analysis output to this directory.
        maxEvents: The maximum number of events per run to analyze.
        
        """
#### ????
