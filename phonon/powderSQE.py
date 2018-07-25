
import numpy as np


def compute(
        omega, pols, positions, Q_basis, qgrid_shape,
        nbranches=3, max_hkl=1,
        Q_bins = np.arange(0, 11, 0.1), E_bins = np.arange(0, 50, 0.5),
):
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
    bins = Q_bins, E_bins
    for tau_hkl in tau1:
        tau_hkl = np.array(tau_hkl)
        print tau_hkl
        for ibr in range(nbranches):
            omega0 = omega[:, :, :, ibr][indexes]
            Q_hkl = tau_hkl + q_hkl[indexes]
            Q_cart = np.dot(Q_hkl, Q_basis)  # nQ, 3
            Q_mag = np.linalg.norm(Q_cart, axis=-1) # nQ
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

