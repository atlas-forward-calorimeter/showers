import os

import context
from analysis.run import Run

# Force parameters to match the desired test case.
Run.tubeZ = 35 / 2
Run.tubeMiddleZ = 8 / 2

# Analyze a test event.
run = Run(
    os.path.join(context.testsDir, 'data', 'run1'),
    outDirectory=False
)
run.go()

# Test total energy deposit.
assert run.meanEdep == 6306.383002036315, run.meanEdep
assert run.eDepSigma == 130.58126171142476, run.eDepSigma
# Test middle tube energy deposit.
assert run.meanMiddleEdep == 78.635272187729, run.meanMiddleEdep
assert run.middleEdepSigma == 5.6851438688610045, run.middleEdepSigma

print('Run tests passed.')
