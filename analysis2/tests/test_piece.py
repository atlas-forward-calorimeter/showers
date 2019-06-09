from context import piece

class TestPiece(piece.Piece):
    def __init__(self):
        super().__init__(isAtom=True)

    def subAnalyses(self):
        raise RuntimeError(
            "Bottom level data (atoms) shouldn't try to loop through"
            "smaller pieces."
        )


# Test that the atoms don't try to loop through smaller pieces.
testPiece = TestPiece()
testPiece.go()


"""
halleuljah <showers.analysis2.run.Run object at 0x0000015DF1BBA588>
ok so rn I <showers.analysis2.multirun.MultiRun object at 0x0000015DF1B8E0F0> am accessing <showers.analysis2.run.Run object at 0x0000015DF1BBA588>
halleuljah <showers.analysis2.run.Run object at 0x0000015DF1BBA5C0>
ok so rn I <showers.analysis2.multirun.MultiRun object at 0x0000015DF1B8E0F0> am accessing <showers.analysis2.run.Run object at 0x0000015DF1BBA5C0>
halleuljah <showers.analysis2.run.Run object at 0x0000015DF1B8EA20>
ok so rn I <showers.analysis2.multirun.MultiRun object at 0x0000015DF1B8E0F0> am accessing <showers.analysis2.run.Run object at 0x0000015DF1B8EA20>
halleuljah <showers.analysis2.run.Run object at 0x0000015DF1BBA5C0>
ok so rn I <showers.analysis2.multirun.MultiRun object at 0x0000015DF1B8E0F0> am accessing <showers.analysis2.run.Run object at 0x0000015DF1BBA5C0>
halleuljah <showers.analysis2.multirun.MultiRun object at 0x0000015DF1B8E0F0>
"""
