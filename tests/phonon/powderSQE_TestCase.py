#!/usr/bin/env python

skip = True # this test needs phonopy

import unittest, os, sys, shutil, numpy as np, histogram as H, histogram.hdf as hh

here = os.path.dirname(__file__)

class TestCase(unittest.TestCase):
    
    def test(self):
        from mcvine.phonon.powderSQE import compute
        from mccomponents.sample.idf import Polarizations, Omega2, units
        # read data
        pols = Polarizations.read('/SNS/users/linjiao/simulations/samples/Al/phonons-grid100/Polarizations')[1]
        pols.shape = 100,100,100, 3, 1, 3, 2
        _pols = pols
        pols = _pols[:, :, :, :, :, :, 0] + 1j * _pols[:, :, :, :, :, :, 1]
        omega2 = Omega2.read('/SNS/users/linjiao/simulations/samples/Al/phonons-grid100/Omega2')[1]
        omega2.shape = 100,100,100, 3
        omega = omega2**.5 * units.hertz2mev
        #
        # from Qgridinfo
        b1=   -1.555243888,   -1.555243888,    1.555243888
        b2=    1.555243888,    1.555243888,    1.555243888
        b3=   -1.555243888,    1.555243888,   -1.555243888
        #
        Q_basis = np.array([b1, b2, b3])
        # from POSCAR
        basis = np.array(
            [[0, 2.02, 2.02],
             [2.02, 0, 2.02],
             [2.02, 2.02, 0]
            ])
        # atom positions
        positions = [
                [0, 0, 0],
            ]
        positions = np.dot(positions, basis)
        #
        Qbb, Ebb, I = compute(omega, pols, positions, Q_basis, nbranches=3, max_hkl=1)
        IQEhist = H.histogram(
            'IQE',
            (H.axis('Q', boundaries=Qbb, unit='1./angstrom'),
             H.axis('E', boundaries=Ebb, unit='meV')),
            data=I)
        hh.dump(IQEhist, 'Al-iqe.h5')
        return


if __name__ == '__main__': unittest.main()
