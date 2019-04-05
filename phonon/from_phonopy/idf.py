# -*- Python -*-
#
# Jiao Lin <jiao.lin@gmail.com>
#

"""
utils to create phonon data in IDF format from phonopy data
"""

from . import call_phonopy
from phonopy.interface import vasp

import numpy as np


def make_all(
        species, supercell_dims=[5,5,5], qgrid_dims=[51,51,51], fix_pols_phase=True,
        force_constants='FORCE_CONSTANTS', poscar='POSCAR',
):
    """ compute all phonon data needed by the single crystal phonon kernel.
    
    Inputs:
      - species: list of atomic species
      - supercell_dims: supercell dimensions, eg "5 5 5". should be consistent with the FORCE_CONSTANTS file
      - qgrid_dims: Q grid dimensions
    
    Output:
      - Omega2
      - Polarizations
      - Qgridinfo
      - DOS
    
    """
    make_omega2_pols(
        species, supercell_dims, qgrid_dims, 
        force_constants=force_constants, poscar=poscar, fix_phase = fix_pols_phase
    )
    make_Qgridinfo(qgrid_dims, species, poscar)
    from ..dos import fromOmaga2
    fromOmaga2()
    xtal_xyz = ''.join(species) + '.xyz'
    make_crystal_xyz(xtal_xyz, species, poscar)
    return


def make_omega2_pols(
        species, supercell_dims=[5,5,5], qgrid_dims=[51,51,51],
        force_constants='FORCE_CONSTANTS', poscar='POSCAR',
        fix_phase=True
):
    """compute polarizations
    
    Inputs
      - force_constants: path of the FORCE_CONSTANTS file
      - poscar: path of the POSCAR file
      - species: list of atomic species
      - supercell_dims: supercell dimensions, eg "5 5 5". should be consistent with the FORCE_CONSTANTS file
      - qgrid_dims: Q grid dimensions
    
    Output:
      - Omega2
      - Polarizations
    """
    print "* Constructing Q array"
    qgrid_dims = np.array(qgrid_dims, dtype=float)
    delta = 1./(qgrid_dims-1)
    Qx = np.arange(0, 1.+delta[0]/2, delta[0])
    Qy = np.arange(0, 1.+delta[1]/2, delta[1])
    Qz = np.arange(0, 1.+delta[2]/2, delta[2])
    Qs = []
    for qx in Qx:
        for qy in Qy:
            for qz in Qz:
                Qs.append([qx,qy,qz])
                continue
    Qs =  np.array(Qs)

    # read force_constants
    from phonopy import file_IO
    force_constants=file_IO.parse_FORCE_CONSTANTS(force_constants)

    # !!! only need one symbol per specie
    # !!! follow vasp convention !!!
    supercell_matrix = np.diag(supercell_dims)
    print "* Calling phonopy to compute eigen values and eigen vectors"
    qvecs, freq, pols = call_phonopy.onGrid(species, Qs, supercell_matrix, force_constants, freq2omega=1, poscar=poscar)
    
    print "* Writing out freqencies"
    from mccomponents.sample.idf import Omega2, Polarizations
    freq[freq<0] = 0
    # min = np.min(freq)
    # if min < 0: freq += -min
    omega2 = freq**2 * 1e24 * (2*np.pi)**2
    Omega2.write(omega2)

    # phase factor for pols
    from phonopy.interface import vasp
    print "* Fixing and writing out polarizations"
    nq, nbr, natoms, three = pols.shape
    assert three is 3
    if fix_phase:
        print "* Fixing phase: exp(i Q.d)"
        atoms = vasp.read_vasp(poscar, species)
        positions = atoms.get_scaled_positions()
        # correct polarization vectors
        # the phase factor is needed. see the notebook tests/phonon/phase-factor.ipynb
        # c.c. is needed at the very end because of another convention difference between phonopy and pybvk codes.
        # without c.c., we have to use exp(-j Q dot r) instead of exp(j Q dot r)
        for iatom in range(natoms):
            qdotr = np.dot(Qs, positions[iatom]) * 2 * np.pi
            phase = np.exp(1j * qdotr)
            pols[:, :, iatom, :] *= phase[:, np.newaxis, np.newaxis]
            norms = np.linalg.norm(pols, axis=-1)
            pols/=norms[:, :, :, np.newaxis]
            continue
        pols = np.conj(pols)
    Polarizations.write(pols)
    return


def make_Qgridinfo(qgrid_dims, species, poscar='POSCAR'):
    """Create Q gridinfo file in IDF format
    """
    from phonopy.interface import vasp
    import numpy as np
    atoms = vasp.read_vasp(poscar, species)
    # atoms.cell is [a1,a2,a3]; inv(atoms.cell) is [b1,b2,b3].T; we want [b1,b2,b3]; hence the .T
    reci_cell = np.linalg.inv(atoms.cell).T * 2*np.pi
    # output
    ostream  = open("Qgridinfo", 'wt')
    # reciprocal cell
    for i in range(3):
        ostream.write("b%d = %s\n" % (i+1, list(reci_cell[i])))
        continue
    # grid
    N = qgrid_dims
    for i in range(3):
        ostream.write("n%d = %s\n" % (i+1, N[i]))
        continue
    ostream.close()
    return


def make_crystal_xyz(outpath, atom_chemical_symbols, poscar='POSCAR'):
    from phonopy.interface import vasp
    atoms = vasp.read_vasp(poscar, atom_chemical_symbols)
    # # of atoms
    lines = [str(len(atoms.get_chemical_symbols()))]
    # lattice
    c = atoms.cell.copy()
    c.shape = -1
    lines.append("\t".join(map(str, c)))
    # atoms and positions
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_scaled_positions()
    for s, p in zip(symbols, positions):
        line = "%s\t%s" % (s, ' '.join(map(str, p)))
        lines.append(line)
        continue
    text = '\n'.join(lines)
    # output stream
    ostream = open(outpath, 'wt')
    ostream.write(text)
    return


# End of file 
