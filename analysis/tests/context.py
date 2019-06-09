"""Give path context to the tests, so they can run from anywhere."""

import os
import sys

testsDir = os.path.abspath(os.path.dirname(__file__))

sys.path.insert(0, os.path.join(testsDir, '..', '..'))
