"""
The main idea is to follow the implementation in powder coherent
single-phonon inelastic kernel.

Converted from https://github.com/mcvine/phonon/blob/0abbc2613a41471fbb51a96f0be9e8a6de89934b/tests/phonon/explore-random-sampling-big-Qspace-box-no-importance-sampleing-along-Qmag.ipynb
"""

import numpy as np, os, sys, glob, histogram as H, histogram.hdf as hh
from mcvine.phonon.from_phonopy import onGrid


def from_FORCE_CONSTANTS(
        datadir,
        Ei = 300., # meV
        T=300., # kelvin
        doshist=None, # DOS histogram
        supercell = (6,6,1),
        Q_bins = np.arange(0, 11, 0.1), E_bins = np.arange(0, 50, 0.5),
        workdir = None, N=int(1e6),
        include_multiphonon=True, scale_multiphonon=1.0, max_det_angle=135.,
        extend_to_negative_E=False,
):
    '''This implementation only uses phonopy. The implementation in IDF uses
    dispersion data stored in IDF format.
    '''
    # create and change to workdir
    if workdir is None:
        import tempfile
        workdir = tempfile.mkdtemp()
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    saved = os.path.abspath('.')
    os.chdir(workdir)
    # material structure
    sc_mat = np.diag(supercell)
    from phonopy import file_IO
    force_constants=file_IO.parse_FORCE_CONSTANTS(os.path.join(datadir, 'FORCE_CONSTANTS'))
    poscar = os.path.join(datadir, 'POSCAR')
    from mcvine.phonon import from_phonopy
    from_phonopy.make_crystal_xyz('structure.xyz', poscar)
    from sampleassembly.crystal.ioutils import xyzfile2unitcell
    uc = xyzfile2unitcell('structure.xyz')
    # reciprocal basis
    # rows are reciprocal lattice vectors
    Q_basis = np.linalg.inv(uc.lattice.base).T*np.pi*2
    np.allclose(np.dot(uc.lattice.base, Q_basis.T), np.eye(3)*2*np.pi)
    # generate Q points
    max_Q = Q_bins[-1]
    Qmag_p3 = np.random.rand(N)*max_Q**3
    Qmag = Qmag_p3**(1./3)
    cos_theta = np.random.rand(N)*2-1  # -1 -- 1
    phi = np.random.rand(N) * 2 * np.pi
    sin_theta = np.sqrt(1-cos_theta*cos_theta)
    sin_phi = np.sin(phi); cos_phi = np.cos(phi)
    Qx = Qmag*sin_theta*cos_phi
    Qy = Qmag*sin_theta*sin_phi
    Qz = Qmag*cos_theta
    Qpoints = np.array([Qx,Qy,Qz]).T
    # in reciprocal units
    hkls = np.dot(Qpoints, uc.lattice.base.T)/(2*np.pi)
   # run phonopy
    sys.stdout.write("Running phonopy...  "); sys.stdout.flush()
    qs, freqs, pols = onGrid(
        hkls, sc_mat, 
        force_constants = force_constants,
        poscar = os.path.join(datadir, 'POSCAR'),
        freq2omega=1
    )
    os.chdir(saved)
    sys.stdout.write("Done"); sys.stdout.flush()
    good = freqs > 0
    # omega2 = freqs**2 * 1e24 * (2*np.pi)**2
    omega = freqs*1e12*2*np.pi
    from mccomponents.sample.idf import units
    omega *= units.hertz2mev
    from ._calc import (calcIQE, apply_corrections, apply_dynamical_range)
    Qbb, Ebb, I = calcIQE(uc, omega, pols, Qpoints, Q_bins, E_bins)
    if extend_to_negative_E:
        from ._calc import extend_to_negative_E
        # now Ebb includes negatives
        Qbb, Ebb, I = extend_to_negative_E(Qbb, Ebb, I, T)
    # additional corrections
    symbols = [a.element for a in uc]
    from ..atomic_scattering import AtomicScattering
    masses = [AtomicScattering(s).mass for s in symbols]; average_mass = np.mean(masses)
    IQEhist = apply_corrections(I, Qbb, Ebb, N, average_mass, uc, doshist, T, Ei)
    IQEhist = apply_dynamical_range(IQEhist, Ei, max_det_angle)
    if include_multiphonon:
        from ._calc import multiphononSQE
        mphIQE = multiphononSQE(
            T=T, doshist=doshist, mass=average_mass, Q_bins=Q_bins, E_bins=Ebb)
        symbols = [a.element for a in uc]
        total_xs = sum(AtomicScattering(s).sigma() for s in symbols)
        return IQEhist, mphIQE * (total_xs/4/np.pi,0)
    return IQEhist

