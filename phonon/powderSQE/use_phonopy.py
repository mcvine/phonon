"""
The main idea is to follow the implementation in powder coherent
single-phonon inelastic kernel.

Converted from https://github.com/mcvine/phonon/blob/0abbc2613a41471fbb51a96f0be9e8a6de89934b/tests/phonon/explore-random-sampling-big-Qspace-box-no-importance-sampleing-along-Qmag.ipynb
"""

import numpy as np, os, sys, glob, histogram as H, histogram.hdf as hh
from multiphonon.sqe import plot as plotsqe, interp as interp_sqe
from multiphonon.forward.phonon import kelvin2mev
import multiphonon
from multiphonon.forward import phonon
from mcvine.phonon.from_phonopy import call_phonopy
from mcni.utils import conversion


def from_FORCE_CONSTANTS(
        datadir,
        Ei = 300., # meV
        T=300., # kelvin
        doshist=None, # DOS histogram
        mass = 12, # hack
        species = ['C'], supercell = (6,6,1),
        Q_bins = np.arange(0, 11, 0.1), E_bins = np.arange(0, 50, 0.5),
        workdir = None, N=int(1e6),
        include_multiphonon=True, scale_multiphonon=1.0, max_det_angle=135.,
):
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
    from_phonopy.make_crystal_xyz('structure.xyz', species, poscar)
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
    qs, freqs, pols = call_phonopy.onGrid(
        species, hkls, sc_mat, 
        force_constants = force_constants,
        poscar = os.path.join(datadir, 'POSCAR'),
        freq2omega=1
    )
    sys.stdout.write("Done"); sys.stdout.flush()
    good = freqs > 0
    # omega2 = freqs**2 * 1e24 * (2*np.pi)**2
    omega = freqs*1e12*2*np.pi
    from mccomponents.sample.idf import units
    omega *= units.hertz2mev
    nq, nbr, natoms, three = pols.shape
    assert three is 3
    # need atom positions
    from phonopy.interface import vasp
    atoms = vasp.read_vasp(poscar, species)
    positions = atoms.get_scaled_positions()
    # normalize polarization vectors
    for iatom in range(natoms):
        qdotr = np.dot(qs, positions[iatom]) * 2 * np.pi
        phase = np.exp(1j * qdotr)
        pols[:, :, iatom, :] *= phase[:, np.newaxis, np.newaxis]
        norms = np.linalg.norm(pols, axis=-1)
        pols/=norms[:, :, :, np.newaxis]
        continue
    # Compute
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
        exp_Q_dot_d = np.exp(-1j * np.dot(good_hkls, positions.T) * 2*np.pi) # nQ, natoms
        pols1 = pols[good1, ibr, :, :] # nQ, natoms, 3   
        Q_dot_pol = np.sum(np.transpose(pols1, (1,0,2)) * Q_cart, axis=-1).T # nQ, natoms
        # 
        F = np.sum(exp_Q_dot_d * Q_dot_pol, axis=-1) # nQ      
        M = np.abs(F)**2 / omega_good # nQ  
        I1, Qbb, Ebb = np.histogram2d(Qmag_good, omega_good, bins=bins, weights=M)
        I = I + I1
        continue
    I *= 1./mass * conversion.k2e(1.)
    Q = (Qbb[1:] + Qbb[:-1])/2
    # correction 1
    ki = conversion.e2k(Ei)
    E = (Ebb[1:] + Ebb[:-1])/2
    Ef = Ei - E
    kf = conversion.e2k(Ef)
    v0 = uc.lattice.volume
    beta = 1./(T*kelvin2mev)
    thermal_factor = 1./2 * (1./np.tanh(E/2*beta) + 1)
    correction_E = thermal_factor * (2*np.pi)**3/v0 / (2*ki*kf)
    correction_Q = 1./Q
    correction = np.outer(correction_Q, correction_E)
    I*=correction
    # # Correction 2 - DW
    dos_e = doshist.E if hasattr(doshist, 'E') else doshist.energy
    dos_dE = dos_e[1]-dos_e[0]
    dos_g = doshist.I / np.sum(doshist.I) /dos_dE
    DW2 = phonon.DWExp(Q, M=mass, E=dos_e,g=dos_g, beta=beta, dE=dos_e[1]-dos_e[0])
    I*=np.exp(-DW2)[:, np.newaxis]
    IQEhist = H.histogram(
        'IQE',
        (H.axis('Q', boundaries=Qbb, unit='1./angstrom'),
         H.axis('E', boundaries=Ebb, unit='meV')),
        data=I)
    ## Normalize
    # first normalize by "event" count and make it density
    dQ = Q[1] - Q[0]
    dE = E[1] - E[0]
    IQEhist.I /= N
    IQEhist.I /= dQ*dE
    # Simulate "Vanadium" data
    def norm_at(E, Qaxis, Ei):
        "compute normalization(Q) array for the given energy transfer"
        ki = conversion.e2k(Ei)
        Ef = Ei - E
        kf = conversion.e2k(Ef)
        Qmin = np.abs(ki - kf)
        Qmax = ki + kf
        Qmiddle = (Qmin + Qmax)/2
        v_middle = 1./(Qmax-Qmin)
        slope = v_middle / Qmiddle
        rt = np.zeros(Qaxis.shape)
        in_range = (Qaxis<Qmax) * (Qaxis>Qmin)
        rt[in_range] = (Qaxis * slope)[in_range]
        return rt
    norm_hist = IQEhist.copy()
    norm_hist.I[:] = 0
    for E_ in norm_hist.E:
        norm_hist[(), E_].I[:] = norm_at(E_, Q, Ei)
    # normalize
    IQEhist = IQEhist/norm_hist
    #
    ## Dynamical range
    DR_Qmin = ki-kf
    DR_Qmax = ((ki*ki + kf*kf - 2*ki*kf*np.cos(max_det_angle*np.pi/180)))**.5
    I = IQEhist.I
    for iE, (E1, Qmin1, Qmax1) in enumerate(zip(E, DR_Qmin, DR_Qmax)):
        I[Q<Qmin1, iE] = np.nan
        I[Q>Qmax1, iE] = np.nan
    #
    if include_multiphonon:
        ## multiphonon
        mpsqe = phonon.sqehist(
            dos_e, dos_g, Qmin=Q[0], Qmax=Q[-1]+dQ/2., dQ=dQ, T=T, M=mass, N=5, starting_order=2, Emax=E[-1]*3)
        mpsqe1 = interp_sqe(mpsqe, E)
        mpsqe1.I *= scale_multiphonon; mpsqe1.E2 *= scale_multiphonon*scale_multiphonon
    else:
        mpsqe1 = (0,0)
    os.chdir(saved)
    # return IQEhist
    return mpsqe1 + IQEhist*(8*np.pi,0)

def multiphononSQE(
        T=300., # kelvin
        doshist=None, # DOS histogram
        mass = 12, # hack
        Q_bins = np.arange(0, 11, 0.1), E_bins = np.arange(0, 50, 0.5),
):
    dos_e = doshist.E if hasattr(doshist, 'E') else doshist.energy
    dos_dE = dos_e[1]-dos_e[0]
    dos_g = doshist.I / np.sum(doshist.I) /dos_dE
    Qbb = Q_bins
    Q = (Qbb[1:] + Qbb[:-1])/2
    dQ = Qbb[1]-Qbb[0]
    Ebb = E_bins
    E = (Ebb[1:] + Ebb[:-1])/2
    mpsqe = phonon.sqehist(
        dos_e, dos_g, Qmin=Q[0], Qmax=Q[-1]+dQ/2., dQ=dQ, T=T, M=mass, N=5, starting_order=2, Emax=E[-1]*3)
    mpsqe1 = interp_sqe(mpsqe, E)
    return mpsqe1
