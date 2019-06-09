"""Give path context to the tests, so they can run from anywhere."""

import os
import sys
sys.path.insert(
    0, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
)

import piece
import multirun
