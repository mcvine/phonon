#!/usr/bin/env python

skip = True # this test needs phonopy

import unittest, os, glob, sys, shutil, numpy as np, histogram as H, histogram.hdf as hh

here = os.path.abspath(os.path.dirname(__file__))

class TestCase(unittest.TestCase):

    def test2(self):
        datadir = os.path.join(here, '..', '..', 'data', 'graphite')
        doshist = hh.load(os.path.join(datadir, 'exp_DOS.h5'))
        from mcvine.phonon.powderSQE.use_phonopy import from_FORCE_CONSTANTS
        IQEhist = from_FORCE_CONSTANTS(
            datadir,
            Ei = 300., # meV
            T=300., # kelvin
            doshist = doshist, # DOS histogram
            mass = 12, # hack
            species = ['C'], supercell = (6,6,1),
            Q_bins = np.arange(0, 23, 0.1), E_bins = np.arange(0, 250, 1.),
            workdir = '_tmp.powderSQE', N=int(1e5)
        )
        hh.dump(IQEhist, 'graphite-Ei_300-T_300.h5')
        return


if __name__ == '__main__': unittest.main()
