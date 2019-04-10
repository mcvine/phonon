#!/usr/bin/env python

skip = True # this test needs phonopy

import unittest, os, glob, sys, shutil, numpy as np, histogram as H, histogram.hdf as hh

here = os.path.abspath(os.path.dirname(__file__))

class TestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.datadir = work = os.path.join(here, '_tmp.IDF-Si-phonons')
        if os.path.exists(work):
            shutil.rmtree(work)
        from mcvine import deployment_info
        from mcvine.resources import sample
        src = sample(name='Si', temperature='100K', shape='dummy')
        src = os.path.join(os.path.dirname(src), 'phonons', 'vasp-phonopy')
        shutil.copytree(src, work)
        from mcvine.phonon.from_phonopy import make_all
        saved = os.path.abspath('.')
        os.chdir(work)
        make_all(
            species=['Si'], supercell_dims=[5,5,5],
            qgrid_dims=[51,51,51],
            fix_pols_phase=True,
            force_constants='FORCE_CONSTANTS', poscar='POSCAR',
        )
        os.chdir(saved)
        return
    
    def test2(self):
        datadir = self.datadir
        from mcvine.phonon.powderSQE.IDF import from_data_dir
        IQEhist = from_data_dir(datadir, max_hkl=10)
        hh.dump(IQEhist, 'Si-iqe-test2.h5')
        return


if __name__ == '__main__': unittest.main()
