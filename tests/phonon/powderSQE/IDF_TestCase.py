#!/usr/bin/env python

skip = True # this test needs phonopy

import unittest, os, glob, sys, shutil, numpy as np, histogram as H, histogram.hdf as hh

here = os.path.abspath(os.path.dirname(__file__))

class TestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.datadir = work = os.path.join(here, '_tmp.powderSQE-Si-phonons')
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
    
    def test(self):
        datadir = self.datadir
        from mcvine.phonon.powderSQE.IDF import compute
        from mccomponents.sample.idf import Polarizations, Omega2, units
        from mcvine.phonon.io import readQgridinfo
        # atoms. read it from xyz files
        xyzfiles = glob.glob(os.path.join(datadir, '*.xyz'))
        if len(xyzfiles)>1: raise NotImplementedError("more than one xyz files")
        xyzf = xyzfiles[0]
        from sampleassembly.crystal.ioutils import xyzfile2unitcell
        uc = xyzfile2unitcell(xyzf)
        nAtoms = len(uc)
        # Qgridinfo
        reci_basis, gridshape = readQgridinfo(os.path.join(datadir, 'Qgridinfo'))
        # read data
        pols = Polarizations.read(os.path.join(datadir, 'Polarizations'))[1]
        ndims = 3
        nbranches = ndims*nAtoms
        pols.shape = gridshape + (nbranches, nAtoms, 3, 2)
        _pols = pols
        pols = _pols[:, :, :, :, :, :, 0] + 1j * _pols[:, :, :, :, :, :, 1]
        omega2 = Omega2.read(os.path.join(datadir, 'Omega2'))[1]
        omega2.shape = gridshape + (nbranches,)
        omega = omega2**.5 * units.hertz2mev
        #
        Q_basis = np.array(reci_basis)
        # from POSCAR
        basis = uc.lattice.base
        # atom positions
        positions = np.array([np.array(a.xyz_cartn) for a in uc])
        #
        Qbb, Ebb, I = compute(omega, pols, positions, Q_basis, gridshape, nbranches=nbranches, max_hkl=10)
        IQEhist = H.histogram(
            'IQE',
            (H.axis('Q', boundaries=Qbb, unit='1./angstrom'),
             H.axis('E', boundaries=Ebb, unit='meV')),
            data=I)
        hh.dump(IQEhist, 'Si-iqe.h5')
        return


    def test2(self):
        datadir = self.datadir
        from mcvine.phonon.powderSQE.IDF import from_data_dir
        IQEhist = from_data_dir(datadir, max_hkl=10)
        hh.dump(IQEhist, 'Si-iqe-test2.h5')
        return


if __name__ == '__main__': unittest.main()
