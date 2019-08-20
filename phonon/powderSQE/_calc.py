import numpy as np

def calcIQE(uc, omega, pols, Qpoints, Q_bins, E_bins):
    """main part of IQE calculation. missing several corrections
    """
    N = len(Qpoints)
    Qmag = np.linalg.norm(Qpoints, axis=-1)
    from ..atomic_scattering import AtomicScattering
    symbols = [a.element for a in uc]
    bs = [AtomicScattering(s).b() for s in symbols]
    masses = [AtomicScattering(s).mass for s in symbols]
    positions = np.array([a.xyz for a in uc])
    bovermass = np.array(bs)/np.sqrt(masses)
    #
    import tqdm
    bins = Q_bins, E_bins
    # in reciprocal units
    hkls = np.dot(Qpoints, uc.lattice.base.T)/(2*np.pi)
    I = 0
    nQ, nbr, natoms, three = pols.shape
    assert three == 3
    good = omega > 0
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
    uc_vol = uc.lattice.volume; ruc_vol = (2*np.pi)**3/uc_vol
    max_Q = Q_bins[-1]
    I*=1./N*4./3*np.pi*max_Q**3/ruc_vol
    return Qbb, Ebb, I


def apply_corrections(I, Qbb, Ebb, N, mass, uc, doshist, T, Ei, max_det_angle):
    """apply remaining corrections to IQE. this needs to be combined with calcIQE"""
    from mcni.utils import conversion
    from multiphonon.forward.phonon import kelvin2mev
    from multiphonon.forward import phonon
    import histogram as H
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
    # the additional 1/4pi comes from powder average. See notes
    correction = 1./4/np.pi*np.outer(correction_Q, correction_E)
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
    IQEhist.I /= dQ*dE
    # Simulate "Vanadium" data
    def norm_at(E, Qaxis, Ei):
        "compute normalization(Q) array for the given energy transfer"
        ki = conversion.e2k(Ei)
        Ef = Ei - E
        kf = conversion.e2k(Ef)
        Qmin = np.abs(ki - kf)
        Qmax = ki + kf
        rt = np.zeros(Qaxis.shape)
        in_range = (Qaxis<Qmax) * (Qaxis>Qmin)
        rt[in_range] = (Qaxis *(2*np.pi/ki/kf))[in_range]
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
    return IQEhist

