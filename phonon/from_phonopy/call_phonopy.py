# -*- Python -*-
#
# Jiao Lin <jiao.lin@gmail.com>
#

"""
utils to create phonon data by using phononpy
"""

import numpy as np
from phonopy.interface import vasp
from phonopy.units import VaspToTHz
from phonopy import Phonopy, file_IO

# some constants
THz2meV=4.1357

def onGrid(
        atom_chemical_symbols, qpoints, supercell_matrix, 
        force_constants,
        freq2omega=THz2meV, poscar='POSCAR',
):
    """use phonopy to compute phonon frequencies and polarizations
    for the given Q points on a grid

    force_constants: instance
    """
    
    # set up Si crystal lattice
    bulk = vasp.read_vasp(poscar, atom_chemical_symbols)
    
    # phonopy phonon instance
    phonon = Phonopy(bulk, supercell_matrix, factor=VaspToTHz)
    phonon.generate_displacements(distance=0.01)
    # symmetry = phonon.get_symmetry()
    
    # report
    # print "Space group:", symmetry.get_international_table()
    # phonon.print_displacements()
    # supercells = phonon.get_supercells_with_displacements()

    # set force constants
    phonon.set_force_constants(force_constants)

    # calc band structure
    # . compute
    phonon.set_qpoints_phonon(
        qpoints, is_eigenvectors=True,
        write_dynamical_matrices=False) # , factor=VaspToTHz)

    # output band structure
    # phonon.write_yaml_qpoints_phonon()

    # . get data
    freq, pols = phonon.get_qpoints_phonon()
    freq = freq * freq2omega
    pols = np.transpose(pols, (0, 2, 1))
    pols.shape = pols.shape[:-1] + (-1, 3)
    # pols: Q, branch, atom, xyz
    return qpoints, freq, pols


# End of file 
