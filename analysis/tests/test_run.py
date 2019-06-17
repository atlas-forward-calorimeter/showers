import unittest
import os
import numpy as np

class TestRun(unittest.TestCase):
    def setUp(self):
        import context
        from analysis.run import Run

        # Force parameters to match the desired test case.
        Run.tubeZ = 35 / 2
        Run.tubeMiddleZ = 8 / 2

        # Analyze a test event.
        self.run = Run(
            os.path.join(context.testsDir, 'data', 'run1'),
            outDirectory=False
        )
        self.run.go()
    
    def test_total_energy_deposit(self):
        self.assertAlmostEqual(self.run.meanEdep, 6306.383002036315)
        self.assertAlmostEqual(self.run.eDepSigma, 130.58126171142476)
    
    def test_middle_tube_energy_depost(self):
        self.assertAlmostEqual(self.run.meanMiddleEdep, 78.635272187729)
        self.assertAlmostEqual(self.run.middleEdepSigma, 5.6851438688610045)
    
    def test_energy_z_histogram(self):
        testSums = np.loadtxt(
            os.path.join('data', 'testcalcs_output', 'run-zhist.csv'))
        # Assert in numpy and use the result.
        self.assertTrue(np.allclose(self.run.zMeanFullSums, testSums))

    def test_energy_z_sigmas(self):
        testSums = np.loadtxt(
            os.path.join('data', 'testcalcs_output', 'run-zhistsigmas.csv'))
        # Assert in numpy and use the result.
        self.assertTrue(np.allclose(self.run.zFullSumsSigmas, testSums))

    def test_energy_xy_histogram(self):
        testSums = np.loadtxt(
            os.path.join('data', 'testcalcs_output', 'run-xyhist.csv'))
        # Assert in numpy and use the result.
        self.assertTrue(np.allclose(self.run.xyMeanMiddleSums, testSums))
