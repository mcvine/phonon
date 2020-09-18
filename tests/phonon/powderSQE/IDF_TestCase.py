#!/usr/bin/env python

skip = True # this test needs phonopy
plot = False

import unittest, os, glob, sys, shutil, numpy as np, histogram as H, histogram.hdf as hh
if plot:
    from multiphonon import sqe as mpsqe
    from matplotlib import pyplot as plt

here = os.path.abspath(os.path.dirname(__file__))

class TestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.datadir = work = os.path.join(here, '_tmp.IDF-graphite-phonons')
        if os.path.exists(work):
            os.chdir(work)
            return
        datasrc = os.path.join(here, '..', '..', 'data', 'graphite')
        shutil.copytree(datasrc, work)
        os.chdir(work)
        from mcvine.phonon.from_phonopy import make_all
        make_all(
            # species=['C'], supercell_dims=[6,6,1],
            qgrid_dims=[51,51,51],
            fix_pols_phase=True,
            force_constants='FORCE_CONSTANTS', poscar='POSCAR', sposcar='SPOSCAR',
        )
        return
    
    def test(self):
        datadir = self.datadir
        doshist = hh.load(os.path.join(datadir, 'exp_DOS.h5'))
        import mcvine.phonon.powderSQE.IDF as psidf
        disp = psidf.disp_from_datadir(datadir)
        IQEhist, mphhist = psidf.from_data_dir(
            datadir=datadir,
            disp=disp,
            N = int(1e6),
            Q_bins=np.arange(0, 23, 0.1), E_bins=np.arange(0,250,1.),
            doshist=doshist,
            T=300., Ei=300., max_det_angle=140.,
            include_multiphonon=True,
        )
        IQEhist += mphhist
        hh.dump(IQEhist, 'graphite-allphonon-Ei_300-T_300-IDF.h5')
        expected = hh.load(os.path.join(here, 'saved_results/graphite-allphonon-Ei_300-T_300-IDF.h5'))
        max = np.nanmax(expected.I)
        reldiff = IQEhist-expected
        reldiff.I/=max; reldiff.E2/=max*max
        Nbigdiff = (np.abs(reldiff.I)>0.03).sum()
        Ngood = (IQEhist.I==IQEhist.I).sum()
        Ntotal = IQEhist.size()
        self.assertTrue(Ngood*1./Ntotal>.65)
        self.assertTrue(Nbigdiff*1./Ngood<.10)
        if plot:
            plt.figure(figsize=(6,3))
            max = np.nanmax(IQEhist.I)
            median = np.nanmedian(IQEhist.I[IQEhist.I>0])
            mpsqe.plot(IQEhist); plt.clim(0, median*3)
            # plt.subplot(1,2,1);  mpsqe.plot(IQEhist); plt.clim(0, max/50)
            # plt.subplot(1,2,2);  mpsqe.plot(expectedIQEhist); plt.clim(0, max/50)
            plt.show()        
        return



if __name__ == '__main__': unittest.main()
