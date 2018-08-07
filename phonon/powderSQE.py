"""
The main idea is to follow the formula in Squires.

The integration is over #\tau$ points. 
For each $\tau$, a grid of q points is considered.
The phonon data (energies and polarizations) were
already calculated for those q points.
"""

import numpy as np, os, glob, histogram as H


def from_phonon_data_dir(
        datadir,
        max_hkl=1,
        Q_bins = np.arange(0, 11, 0.1), E_bins = np.arange(0, 50, 0.5),
):
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
    # read data and reformat
    pols = Polarizations.read(os.path.join(datadir, 'Polarizations'))[1]
    ndims = 3
    nbranches = ndims*nAtoms
    pols.shape = gridshape + (nbranches, nAtoms, 3, 2)
    _pols = pols
    pols = _pols[:, :, :, :, :, :, 0] + 1j * _pols[:, :, :, :, :, :, 1]
    # normalize pol vectors
    norms = np.linalg.norm(pols, axis=-1)
    pols/=norms[:, :, :, :, :, np.newaxis]
    # read omegas
    omega2 = Omega2.read(os.path.join(datadir, 'Omega2'))[1]
    omega2.shape = gridshape + (nbranches,)
    omega = omega2**.5 * units.hertz2mev
    # Q basis
    Q_basis = np.array(reci_basis)
    # real space basis
    basis = uc.lattice.base
    # atom positions
    positions = np.array([np.array(a.xyz_cartn) for a in uc])
    #
    Qbb, Ebb, I = compute(
        omega, pols, positions, Q_basis, gridshape,
        nbranches=nbranches, max_hkl=max_hkl,
        Q_bins = Q_bins, E_bins = E_bins,
    )
    IQEhist = H.histogram(
        'IQE',
        (H.axis('Q', boundaries=Qbb, unit='1./angstrom'),
         H.axis('E', boundaries=Ebb, unit='meV')),
        data=I)
    return IQEhist

    
def compute(
        omega, pols, positions, Q_basis, qgrid_shape,
        nbranches=3, max_hkl=1,
        Q_bins = np.arange(0, 11, 0.1), E_bins = np.arange(0, 50, 0.5),
):
    """
    omega: shape=(Nq1,Nq2,Nq2,Nbr)
    pols: complex. shape=(Nq1,Nq2,Nq2,Nbr, nAtoms,3)
    positions: shape=(nAtoms, 3)
    Q_basis: shape=(3, 3)
    qgrid_shape: (Nq1, Nq2, Nq3)
    """
    I = 0
    Nq1, Nq2, Nq3 = qgrid_shape
    q_hkl = np.mgrid[0:1.00001:1./(Nq1-1), 0:1.00001:1./(Nq2-1), 0:1.00001:1./(Nq3-1)]
    q_hkl = np.transpose(q_hkl, (1,2,3,0))
    # print q_hkl.shape
    # create indexes array
    ijk = np.mgrid[:Nq1-1, :Nq2-1, :Nq3-1]
    ijk.shape = 3, -1
    ijk = ijk.T
    indexes = tuple(ijk.T)
    # tau array
    tau1 = np.mgrid[-max_hkl:max_hkl+1, -max_hkl:max_hkl+1, -max_hkl:max_hkl+1]
    tau1 = tau1.copy()
    tau1.shape = 3, -1
    tau1 = tau1.T
    # max Q of requested Q axis
    max_Q = Q_bins[-1]
    # compute max Q for the first reciprocal unit cell
    max_Q_1uc = np.max(np.linalg.norm(np.dot(q_hkl[indexes], Q_basis), axis=-1))
    #
    bins = Q_bins, E_bins
    for tau_hkl in tau1:
        tau_hkl = np.array(tau_hkl)
        # print tau_hkl,
        # no need to include this tau point if it is too far away from origin
        Qtau_cart = np.dot(tau_hkl, Q_basis)
        Qtau_mag = np.linalg.norm(Qtau_cart)
        if Qtau_mag - max_Q_1uc > max_Q:
            # print "skipped"
            continue
        else:
            print tau_hkl
        for ibr in range(nbranches):
            Q_hkl = tau_hkl + q_hkl[indexes]
            Q_cart = np.dot(Q_hkl, Q_basis)  # nQ, 3
            Q_mag = np.linalg.norm(Q_cart, axis=-1) # nQ
            #
            omega0 = omega[:, :, :, ibr][indexes]
            #
            exp_Q_dot_d = np.exp(1j * np.dot(Q_cart, positions.T)) # nQ, natoms
            pols1 = pols[:, :, :, ibr, :, :][indexes] # nQ, natoms, 3
            Q_dot_pol = np.sum(np.transpose(pols1, (1,0,2)) * Q_cart, axis=-1).T # nQ, natoms
            #
            F = np.sum(exp_Q_dot_d * Q_dot_pol, axis=-1) # nQ
            M = np.abs(F)**2 # nQ
            I1, Qbb, Ebb = np.histogram2d(Q_mag, omega0, bins=bins, weights=M)
            I = I + I1
        continue
    return Qbb, Ebb, I

