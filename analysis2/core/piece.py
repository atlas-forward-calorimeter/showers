"""The `Piece` base class (a "base base" class)."""

class Piece:
    """A "piece" of analysis.

    An analysis procedure can be thought of as being made up of pieces.
    Each piece can be made up of smaller pieces and make up bigger ones.
    For example, a collection of runs might be one big piece, which is
    made up of runs that are smaller pieces, which are each made up of
    events, the smallest piece.

    """
    def __init__(self, name=None, isAtom=False, parent=None):
        """Initialize everything before looping through smaller pieces.

          `name`: Just to have a nicer sounding name to go by.
        `isAtom`: `True` means this piece is not made of smaller pieces;
                  it is a fundamental unit of analysis.
                  When a piece is an atom, the `go` method will call
                  `analyze` on it but won't try to loop through its
                  smaller pieces.
        `parent`: The next larger piece, if there is one. If not, this
                  piece will be a "top level" piece.
        
        Note: Concrete derived classes will usually want to take some
        input and output parameters in their constructor. For example,
        an input might be a path or URL leading to some data, and an
        output could be similar.
        """
        self.name = name
        self.isAtom = isAtom
        self.parent = parent

    def go(self):
        """Recursively analyze! 
        
        Analyzes each of the next smallest pieces before analyzing
        this piece itself. For example, call `go` on the biggest piece 
        to analyze everything in it.
        
        """
        if not self.isAtom:
            # Loop over next smallest level.
            for piece in self.smallerPieces():
                piece.go()
                # Give this piece access to the smaller one.
                self.analyzeSmaller(piece)  
        # Analyze this level.
        self.analyze()

    def analyze(self):
        """Carry out analysis on the level of this piece.
        
        If this piece is made of smaller ones, `go` will have the 
        smaller pieces do their analysis and call `analyzeSmaller` 
        each time on this piece, then this method will finally 
        be called.
        
        """
    
    def smallerPieces(self):
        """Should be implemented if this piece is not an "atom."
        
        return: An iterator over `Piece` objects of the next lowest 
        level of analysis.

        `self` should be assigned to the parent of the smaller piece in 
        this method, if they will need it.

        """
    
    def analyzeSmaller(self, smallerPiece):
        """A place to carry out analysis, with access to the next
        smallest pieces.

        Called in `go`.
        """
    
    def __repr__(self):
        # TODO: Finish this.
        if self.name:
            return f'<{self.name}>'
        else:
            return super().__repr__()
