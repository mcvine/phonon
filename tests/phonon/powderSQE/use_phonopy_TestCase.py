#!/usr/bin/env python

skip = True # this test needs phonopy
plot = True

import unittest, os, glob, sys, shutil, numpy as np, histogram as H, histogram.hdf as hh
if plot:
    from multiphonon import sqe as mpsqe
    from matplotlib import pyplot as plt

here = os.path.abspath(os.path.dirname(__file__))

class TestCase(unittest.TestCase):

    def test1(self):
        work = os.path.join(here, '_tmp.use_phonopy-Si')
        if os.path.exists(work):
            shutil.rmtree(work)
        from mcvine import deployment_info
        from mcvine.resources import sample
        src = sample(name='Si', temperature='100K', shape='dummy')
        src = os.path.join(os.path.dirname(src), 'phonons', 'vasp-phonopy')
        shutil.copytree(src, work)
        doshist = hh.load(os.path.join(work, 'dos.h5'))
        from mcvine.phonon.powderSQE.use_phonopy import from_FORCE_CONSTANTS
        IQEhist = from_FORCE_CONSTANTS(
            work,
            Ei = 100., # meV
            T=100., # kelvin
            doshist = doshist, # DOS histogram
            mass = 28., # hack
            species = ['Si'], supercell = (5,5,5),
            Q_bins = np.arange(0, 13, 0.1), E_bins = np.arange(0, 95, 0.5),
            workdir = work, N=int(1e5)
        )
        hh.dump(IQEhist, 'silicon-Ei_100-T_100.h5')
        return


    def test2a(self):
        datadir = os.path.join(here, '..', '..', 'data', 'graphite')
        doshist = hh.load(os.path.join(datadir, 'exp_DOS.h5'))
        from mcvine.phonon.powderSQE.use_phonopy import from_FORCE_CONSTANTS
        if plot: N = int(1e5)
        else: N = int(1e6)
        IQEhist = from_FORCE_CONSTANTS(
            datadir,
            Ei = 30., # meV
            T=300., # kelvin
            doshist = doshist, # DOS histogram
            mass = 12, # hack
            species = ['C'], supercell = (6,6,1),
            Q_bins = np.arange(0, 4, 0.04), E_bins = np.arange(0, 30, .2),
            workdir = '_tmp.test2a', N=N, include_multiphonon=False,
            max_det_angle=60.,
        )
        hh.dump(IQEhist, 'graphite-single-phonon-Ei_30-T_300.h5') # save for inspection
        expectedIQEhist = hh.load('saved_results/graphite-single-phonon-Ei_30-T_300-N_3e6.h5')
        expected = expectedIQEhist.I
        # scale it to sth that is easy to get "errorbar"
        N = 3e6
        scale = N/np.nansum(expected)
        expected *= scale; max = np.nanmax(expected)
        this = IQEhist.I; this*=scale
        if plot:
            plt.figure(figsize=(6,3))
            plt.subplot(1,2,1);  mpsqe.plot(IQEhist); plt.clim(0, max/50)
            plt.subplot(1,2,2);  mpsqe.plot(expectedIQEhist); plt.clim(0, max/50)
            plt.show()
        return


    def test2b(self):
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
            Q_bins = np.arange(0, 23, 0.1), E_bins = np.arange(0, 250, 1),
            workdir = '_tmp.powderSQE', N=int(1e5), include_multiphonon=False
        )
        # hh.dump(IQEhist, 'graphite-singlephonon-Ei_300-T_300.h5')
        expected = hh.load('saved_results/graphite-single-phonon-Ei_300-T_300-N_1e7.h5')
        max = np.nanmax(expected.I)
        reldiff = IQEhist-expected
        reldiff.I/=max; reldiff.E2/=max*max
        Nbigdiff = (np.abs(reldiff.I)>0.03).sum()
        Ngood = (IQEhist.I==IQEhist.I).sum()
        Ntotal = IQEhist.size()
        self.assert_(Ngood*1./Ntotal>.65)
        self.assert_(Nbigdiff*1./Ngood<.05)
        return


    def test3(self):
        datadir = os.path.join(here, '..', '..', 'data', 'graphite')
        doshist = hh.load(os.path.join(datadir, 'exp_DOS.h5'))
        from mcvine.phonon.powderSQE.use_phonopy import multiphononSQE
        IQEhist = multiphononSQE(
            T=300., # kelvin
            doshist = doshist, # DOS histogram
            mass = 12, # hack
            Q_bins = np.arange(0, 23, 0.1), E_bins = np.arange(0, 250, 1),
        )
        hh.dump(IQEhist, 'graphite-multiephonon-T_300.h5')
        return


if __name__ == '__main__': unittest.main()
