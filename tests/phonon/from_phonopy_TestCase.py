#!/usr/bin/env python

skip = True # this test needs phonopy

import unittest, os, sys, shutil, numpy as np

here = os.path.dirname(__file__)

class TestCase(unittest.TestCase):
    
    def test(self):
        "mcvine.phonon.from_phonopy: make_IDF_datadir"
        work = '_tmp'
        if os.path.exists(work):
            shutil.rmtree(work)
        from mcvine import deployment_info
        from mcvine.resources import sample
        src = sample(name='Si', temperature='100K', shape='dummy')
        src = os.path.join(os.path.dirname(src), 'phonons', 'vasp-phonopy')
        shutil.copytree(src, work)
        from mcvine.phonon.from_phonopy import make_IDF_datadir
        saved = os.path.abspath('.')
        os.chdir(work)
        make_IDF_datadir(
            supercell_dims=[5,5,5],
            qgrid_dims=[11,11,11],
            fix_pols_phase=True,
            force_constants='FORCE_CONSTANTS', poscar='POSCAR',
        )
        os.chdir(saved)
        return


    def test_make_Qgridinfo(self):
        "mcvine.phonon.from_phonopy: make_Qgridinfo"
        work = '_tmp.make_Qgridinfo'
        poscar=os.path.abspath(os.path.join(here, "../data/graphite/POSCAR"))
        if os.path.exists(work):
            shutil.rmtree(work)
        os.makedirs(work)
        saved = os.path.abspath('.')
        os.chdir(work)
        from mcvine.phonon.from_phonopy import make_Qgridinfo
        make_Qgridinfo((51,51,51), poscar=poscar)
        os.chdir(saved)
        # expected result
        expected = dict(
            b1 = [2.5535175596113087, 1.4742740501120735, 0.0],
            b2 = [0.0, 2.948548100224147, 0.0],
            b3 = [0.0, 0.0, 0.9370895312721231],
            n1 = 51,
            n2 = 51,
            n3 = 51,
            )
        # read result
        res = dict()
        exec(open(os.path.join(work, 'Qgridinfo')).read(), res)
        # compare
        for k,v in expected.iteritems():
            assert np.allclose(v,res[k]), "%s vs %s" % (v, res[k])
        return


if __name__ == '__main__': unittest.main()
