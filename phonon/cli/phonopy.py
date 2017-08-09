#!/usr/bin/env python


from . import mcvine, click

@mcvine.group(help='Commands to run phonopy to compute phonon data in IDF format')
def phonopy():
    return


@phonopy.command()
@click.option("--force-constants", default='FORCE_CONSTANTS', help='path of the FORCE_CONSTANTS file')
@click.option("--poscar", default='POSCAR', help='path of the POSCAR file')
@click.option("--species", default="Si", help='comma-separated list of atomic species')
@click.option("--supercell-dims", default=[5,5,5], help='supercell dimensions, eg "5 5 5". should be consistent with the FORCE_CONSTANTS file', type=int, nargs=3)
@click.option("--qgrid-dims", default=[51, 51, 51], help='Q grid dimensions, eg "51 51 51"', type=int, nargs=3)
def griddisp(force_constants, poscar, species, supercell_dims, qgrid_dims):
    import pdb; pdb.set_trace()
    species = species.split(',')
    from ..from_phonopy import make_all
    make_all(
        species, supercell_dims, qgrid_dims,
        force_constants=force_constants, poscar=poscar, fix_pols_phase=True)
    return
