#!/usr/bin/env python

import unittest
import os
import subprocess
import shlex


class TestGeauxDock(unittest.TestCase):
    def test_a(self):
        os.chdir(os.path.dirname(__file__))
        cmd = '''
./rundock.py ./dock ./1a07C.pdb ./1a07C1.sdf ./1a07C1-0.8.ff
'''

        subprocess.call(shlex.split(cmd))


if __name__ == '__main__':
    unittest.main()
