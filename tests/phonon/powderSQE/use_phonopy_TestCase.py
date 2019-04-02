#!/usr/bin/env python

skip = True # this test needs phonopy

import unittest, os, glob, sys, shutil, numpy as np, histogram as H, histogram.hdf as hh

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
            Q_bins = np.arange(0, 23, 0.1), E_bins = np.arange(0, 250, 1),
            workdir = '_tmp.powderSQE', N=int(1e5), include_multiphonon=False
        )
        hh.dump(IQEhist, 'graphite-singlephonon-Ei_300-T_300.h5')
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
