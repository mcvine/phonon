#!/usr/bin/env python

skip = True # this test needs phonopy

import unittest, os, sys, shutil, numpy as np


class TestCase(unittest.TestCase):
    
    def test(self):
        "mcvine workflow phonon.from_phonopy"
        work = '_tmp'
        if os.path.exists(work):
            shutil.rmtree(work)
        from mcvine import deployment_info
        from mcvine.resources import sample
        src = sample(name='Si', temperature='100K', shape='dummy')
        src = os.path.join(os.path.dirname(src), 'phonons', 'vasp-phonopy')
        shutil.copytree(src, work)
        from mcvine.phonon.from_phonopy import make_all
        os.chdir(work)
        make_all(
            species=['Si'], supercell_dims=[5,5,5], qgrid_dims=[11,11,11], fix_pols_phase=True,
            force_constants='FORCE_CONSTANTS', poscar='POSCAR',
        )
        return


if __name__ == '__main__': unittest.main()
