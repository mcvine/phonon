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
        IQEhist = psidf.from_data_dir(
            datadir=datadir,
            disp=disp,
            N = int(1e5),
            Q_bins=np.arange(0, 4, 0.04), E_bins=np.arange(0,30,.2),
            mass=12., 
            doshist=doshist,
            T=300., Ei=30., max_det_angle=60.,
            include_multiphonon=False,
        )
        # hh.dump(IQEhist, 'graphite-singlephonon-Ei_30-T_300-IDF.h5')
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
