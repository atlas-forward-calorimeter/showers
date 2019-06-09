import os

import context
from analysis.event import Event

# Force parameters to match the desired test case.
Event.tubeZ = 35 / 2
Event.tubeMiddleZ = 8 / 2

# Analyze a test event.
event = Event(
    os.path.join(context.testsDir, 'data', 'run1', 'hits-0.csv'),
    outDirectory=False
)
event.go()

# Test total energy deposit.
assert event.fullEdep == 6175.80174032489, event.fullEdep
# Test middle tube energy deposit.
assert event.middleEdep == 72.950128318868, event.middleEdep

print('Event tests passed.')
