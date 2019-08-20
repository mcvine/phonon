"""
The main idea is to follow the formula in Squires.

The integration is over #\tau$ points. 
For each $\tau$, a grid of q points is considered.
The phonon data (energies and polarizations) were
already calculated for those q points.

This module exists for backward compatibility only.
It computes powder SQE spectrum using phonon data already
calculated in DANSE IDF data format.
For new materials, use "use_phonopy" module, which
requires FORCE_CONSTANTS file.

This module is not very generic: should have generated Q points
randomly and then convert to hkls like what is done
in "use_phonopy" module.
"""

import numpy as np, os, glob, histogram as H

def disp_from_datadir(datadir):
    from mccomponents.sample.phonon import periodicdispersion_fromidf
    dispersion = periodicdispersion_fromidf( datadir )
    from mccomponents.sample import scattererEngine
    return scattererEngine( dispersion )

def from_data_dir(
        datadir=None,
        disp = None,
        N = int(1e6),
        Q_bins = np.arange(0, 11, 0.1), E_bins = np.arange(0, 50, 0.5),
        doshist=None,
        T=300., Ei=100., max_det_angle=135.,
        include_multiphonon=True,
):
    if disp is None:
        from mccomponents.sample.phonon import periodicdispersion_fromidf
        dispersion = periodicdispersion_fromidf( datadir )
        from mccomponents.sample import scattererEngine
        disp = scattererEngine( dispersion )
    #
    poscar = os.path.join(datadir, 'POSCAR')
    from mcvine.phonon import from_phonopy
    from_phonopy.make_crystal_xyz('structure.xyz', poscar)
    from sampleassembly.crystal.ioutils import xyzfile2unitcell
    uc = xyzfile2unitcell('structure.xyz')
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
    #
    omega, pols = _getEsAndPols(disp, Qpoints)
    nQ, nbr, natoms, three = pols.shape
    assert three == 3
    good = omega > 0
    # need atom positions
    from phonopy.interface import vasp
    atoms = vasp.read_vasp(poscar)
    positions = atoms.get_scaled_positions()
    masses = atoms.get_masses()
    average_mass = np.mean(masses)
    symbols = atoms.get_chemical_symbols()
    from ..atomic_scattering import AtomicScattering
    bs = [AtomicScattering(s).b() for s in symbols]
    total_xs = sum(AtomicScattering(s).sigma() for s in symbols)    
    bovermass = np.array(bs)/np.sqrt(masses)
    #
    import tqdm
    bins = Q_bins, E_bins
    I = 0
    for ibr in tqdm.tqdm(range(nbr)):
        good1 = good[:, ibr]
        Qmag_good = Qmag[good1]
        omega_good = omega[good1, ibr]
        Q_cart = Qpoints[good1, :]
        good_hkls = hkls[good1, :]
        #
        exp_Q_dot_d = np.exp(1j * np.dot(good_hkls, positions.T) * 2*np.pi) # nQ, natoms
        pols1 = pols[good1, ibr, :, :] # nQ, natoms, 3
        Q_dot_pol = np.sum(np.transpose(pols1, (1,0,2)) * Q_cart, axis=-1).T # nQ, natoms
        # 
        F = np.sum(exp_Q_dot_d * Q_dot_pol*bovermass, axis=-1) # nQ      
        M = np.abs(F)**2 / omega_good # nQ
        I1, Qbb, Ebb = np.histogram2d(Qmag_good, omega_good, bins=bins, weights=M)
        I = I + I1
        continue
    from mcni.utils import conversion
    I *= conversion.k2e(1.) # canceling unit: hbar^2Q^2/2M/(hbar*omega)
    # M is |b_d/sqrt(M_d)exp(iQ.d)(Q.e)|^2/omega
    # At this moment M (and I) is using b in fm. change it to barn
    I/=100

    # Need to consider the density of the Q points and the volume of the Q space sampled.
    uc_vol = atoms.get_volume(); ruc_vol = (2*np.pi)**3/uc_vol
    I*=1./N*4./3*np.pi*max_Q**3/ruc_vol
    # additional corrections
    from .use_phonopy import _apply_corrections
    IQEhist = _apply_corrections(I, Qbb, Ebb, N, average_mass, uc, doshist, T, Ei, max_det_angle)
    if include_multiphonon:
        from .use_phonopy import multiphononSQE
        mphIQE = multiphononSQE(T=T, doshist=doshist, mass=average_mass, Q_bins=Q_bins, E_bins=E_bins)
        return IQEhist, mphIQE * (total_xs/4/np.pi,0)
    return IQEhist


from mcni.mcnibp import Vector3_double as v3d

def _getEsAndPols(disp, Qs):
    # this only works for mcvine>=1.3.3
    Qs = np.array(Qs, order='C', dtype='double') # C order and double required to work with binding of histogram NdArray
    nQ = len(Qs)
    nbr = disp.nBranches()
    natoms = disp.nAtoms()
    Es = np.zeros((nQ, nbr), dtype='double')
    from mccomponents.sample.phonon.bindings import default
    binding = default()
    disp.energy_arr(binding.ndarray(Qs), binding.ndarray(Es))
    realpol = np.zeros((nQ, nbr, natoms, 3), dtype='double')
    imagpol = np.zeros((nQ, nbr, natoms, 3), dtype='double')
    disp.polarization_arr(binding.ndarray(Qs), binding.ndarray(realpol), binding.ndarray(imagpol))
    pols = realpol + 1j*imagpol
    return Es, pols

def _getEsAndPols_slow(disp, Qs):
    nQ = len(Qs)
    nbr = disp.nBranches()
    natoms = disp.nAtoms()
    Es = np.zeros((nQ, nbr), dtype=float)
    pols = np.zeros((nQ, nbr, natoms, 3), dtype='complex')
    for iQ in range(nQ):
        Q1 = Qs[iQ]
        Q = v3d(float(Q1[0]), float(Q1[1]), float(Q1[2]))
        for ibr in range(nbr):
            Es[iQ, ibr] = disp.energy(ibr, Q)
            for iatom in range(natoms):
                pols[iQ, ibr, iatom, :] = disp.polarization(ibr, iatom, Q)
            continue
        continue
    return Es, pols

