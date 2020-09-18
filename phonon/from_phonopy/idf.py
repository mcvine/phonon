# -*- Python -*-
#
# Jiao Lin <jiao.lin@gmail.com>
#

"""
utils to create phonon data in IDF format from phonopy data
"""

from phonopy.interface import vasp

import numpy as np


def make_all(
        species = None,
        supercell_dims=None, qgrid_dims=[51,51,51], fix_pols_phase=True,
        force_constants='FORCE_CONSTANTS', poscar='POSCAR', sposcar='SPOSCAR',
):
    """ compute all phonon data needed by the single crystal phonon kernel.
    
    Inputs:
      - poscar and sposcar: POSCAR and SPOSCAR file paths. 
      - species: list of atomic species. if None, will infer from POSCAR
      - supercell_dims: supercell dimensions, eg "5 5 5". should be consistent with the FORCE_CONSTANTS file
        if None, will infer from SPOSCAR
      - qgrid_dims: Q grid dimensions
    
    Output:
      - Omega2
      - Polarizations
      - Qgridinfo
      - DOS
    
    """
    # TODO: implement new arg "positions". Right now POSCAR is required to get positions
    # TODO: `make_crystal_xyz` needs revision too
    if species is None:
        # read species from POSCAR
        atoms = vasp.read_vasp(poscar)
        species = set(atoms.get_chemical_symbols())
    else:
        raise NotImplementedError
    if supercell_dims is None:
        supercell_dims = _calc_sc_dims(poscar, sposcar)
    make_omega2_pols(
        supercell_dims, qgrid_dims, 
        force_constants=force_constants, poscar=poscar, fix_phase = fix_pols_phase
    )
    make_Qgridinfo(qgrid_dims, poscar)
    from ..dos import fromOmaga2
    fromOmaga2()
    xtal_xyz = ''.join(species) + '.xyz'
    # species and positions
    make_crystal_xyz(xtal_xyz, poscar)
    return


def _calc_sc_dims(poscar, sposcar):
    cell = vasp.read_vasp(poscar).cell
    sc = vasp.read_vasp(sposcar).cell
    return [_calc_int_scale_factor(v1, v2) for v1,v2 in zip(cell, sc)]

def _calc_scale_factor(v1, v2):
    "v1 and v2 has to be parallel. compute the number s that satisfy v2 = s * v1"
    norm1 = np.linalg.norm(v1); norm2 = np.linalg.norm(v2)
    nv1 = v1/norm1; nv2 = v2/norm2
    if np.linalg.norm(np.cross(nv1, nv2))>1e-3: raise RuntimeError("%s and %s are not parallel" % (v1, v2))
    return norm2/norm1

def _calc_int_scale_factor(v1, v2):
    "v1 and v2 has to be parallel. compute the integer s that satisfy v2 = s * v1"
    s = _calc_scale_factor(v1, v2)
    i = int(s)
    if abs(i-s) > 1e-3:
        raise RuntimeError("(%s)/(%s) is not a whole number" % (v2, v1))
    return i

def make_omega2_pols(
        supercell_dims=[5,5,5], qgrid_dims=[51,51,51],
        force_constants='FORCE_CONSTANTS', poscar='POSCAR',
        fix_phase=True
):
    """compute polarizations
    
    Inputs
      - force_constants: path of the FORCE_CONSTANTS file
      - poscar: path of the POSCAR file
      - supercell_dims: supercell dimensions, eg "5 5 5". should be consistent with the FORCE_CONSTANTS file
      - qgrid_dims: Q grid dimensions
    
    Output:
      - Omega2
      - Polarizations
    """
    print("* Constructing Q array")
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
    print("* Calling phonopy to compute eigen values and eigen vectors")
    from . import onGrid
    qvecs, freq, pols = onGrid(Qs, supercell_matrix, force_constants, freq2omega=1, poscar=poscar)
    
    print("* Writing out freqencies")
    from mccomponents.sample.idf import Omega2, Polarizations
    freq[freq<0] = 0
    # min = np.min(freq)
    # if min < 0: freq += -min
    omega2 = freq**2 * 1e24 * (2*np.pi)**2
    Omega2.write(omega2)
    print("* Writing out polarizations")
    Polarizations.write(pols)
    return


def make_Qgridinfo(qgrid_dims, poscar='POSCAR'):
    """Create Q gridinfo file in IDF format
    """
    from phonopy.interface import vasp
    import numpy as np
    atoms = vasp.read_vasp(poscar)
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


def make_crystal_xyz(outpath, poscar='POSCAR'):
    from phonopy.interface import vasp
    atoms = vasp.read_vasp(poscar)
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
