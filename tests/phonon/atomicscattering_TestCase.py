#!/usr/bin/env python

skip = True # this test needs phonopy

import unittest, os, sys, shutil, numpy as np
here = os.path.dirname(__file__)

from mcvine.phonon.atomic_scattering import AtomicScattering

class TestCase(unittest.TestCase):
    
    def test(self):
        "mcvine.phonon.atomic_scattering"
        H = AtomicScattering("H")
        np.isclose(H.b(), -3.7409)
        H1 = AtomicScattering("H", 1)
        np.isclose(H.b(), -3.7423)
        H2 = AtomicScattering("H", 2)
        np.isclose(H.b(), 6.674)
        return


if __name__ == '__main__': unittest.main()
