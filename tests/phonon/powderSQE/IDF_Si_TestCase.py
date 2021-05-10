#!/usr/bin/env python

skip = True # this test needs phonopy

import unittest, os, glob, sys, shutil, numpy as np, histogram as H, histogram.hdf as hh

here = os.path.abspath(os.path.dirname(__file__))

class TestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.datadir = work = os.path.join(here, '_tmp.IDF-Si-phonons')
        if os.path.exists(work):
            return
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
            supercell_dims=[5,5,5],
            qgrid_dims=[51,51,51],
            fix_pols_phase=True,
            force_constants='FORCE_CONSTANTS', poscar='POSCAR',
        )
        os.chdir(saved)
        return
    
    def test2(self):
        datadir = self.datadir
        from mcvine.phonon.powderSQE.IDF import from_data_dir
        import mcvine.phonon.powderSQE.IDF as psidf
        from mccomponents.sample import phonon as mcphonon
        doshist = mcphonon.read_dos.dos_fromidf(os.path.join(datadir, 'DOS')).doshist
        disp = psidf.disp_from_datadir(datadir)
        IQEhist, mphhist = psidf.from_data_dir(
            datadir=datadir,
            disp=disp,
            N = int(1e6),
            Q_bins=np.arange(0, 14, 0.1), E_bins=np.arange(0,90,.5),
            doshist=doshist,
            T=300., Ei=120., max_det_angle=140.,
            include_multiphonon=True, extend_to_negative_E=True
        )
        IQEhist = IQEhist + mphhist
        hh.dump(IQEhist, 'Si-iqe-test2.h5')
        expected = hh.load(os.path.join(here, 'saved_results/Si-all-phonon-Ei_120-T_300-N_1e6.h5'))
        max = np.nanmax(expected.I)
        reldiff = IQEhist-expected
        reldiff.I/=max; reldiff.E2/=max*max
        Nbigdiff = (np.abs(reldiff.I)>0.03).sum()
        Ngood = (IQEhist.I==IQEhist.I).sum()
        Ntotal = IQEhist.size()
        self.assertTrue(Ngood*1./Ntotal>.65)
        self.assertTrue(Nbigdiff*1./Ngood<.10)
        return

if __name__ == '__main__': unittest.main()
