import unittest
import os
import numpy as np

class TestEvent(unittest.TestCase):
    """Test energy deposit sums and energy-z histograms for events."""
    def setUp(self):
        import context
        from analysis.event import Event

        # Force parameters to match the desired test case.
        Event.tubeZ = 35 / 2
        Event.tubeMiddleZ = 8 / 2

        # Create and analyze a test event.
        self.event = Event(
            os.path.join(context.testsDir, 'data', 'run1', 'hits-0.csv'),
            out_dir=False
        )
        self.event.go()
    
    def test_total_energy_deposit(self):
        self.assertAlmostEqual(self.event.fullEdep, 6175.80174032489)
    
    def test_middle_tube_energy_deposit(self):
        self.assertAlmostEqual(self.event.middleEdep, 72.950128318868)
    
    def test_energy_z_histogram(self):
        testSums = np.loadtxt(
            os.path.join('data', 'testcalcs_output', 'event-hist.csv'))
        # Assert in numpy and use the result.
        self.assertTrue(np.allclose(self.event.zFullSums, testSums))
