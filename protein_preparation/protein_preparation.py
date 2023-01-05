from pymol import cmd
from openbabel import pybel

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdMolAlign

from pdbfixer import PDBFixer
from openmm.app import PDBFile

import MDAnalysis as mda
from MDAnalysis.coordinates import PDB

import random, math

import numpy as np


def separate_protein(code: str):
    """
    Separate protein and ligand from a PDB file.
    :param code: PDB code of the protein.
    :return: None.
    """
    # pdb1 is selected to work only with biological assembly.
    cmd.fetch(code=code, type='pdb')
    cmd.h_add()
    cmd.select(name='receptor', selection='polymer.protein')
    cmd.select(name='ligand', selection='organic')
    cmd.save(filename=code + '_clean.pdb', format='pdb', selection='receptor')
    cmd.save(filename=code + '_lig.mol2', format='mol2', selection='ligand')
    cmd.delete('all')


def get_box(ligand_file: str, extending: float = 6.0, verbose: bool = False):
    """
    Returns PBC box for given ligand and extending.
    :param ligand_file: ligand file name with extension.
    :param extending: extending of the box in Angstrom.
    :param verbose: if True, prints box information with descriptive text else
    returns box in ndarray format suitable for deepchem
    :return: centroid and size of the box.
    """

    cmd.load(filename=ligand_file, format=ligand_file.split('.')[-1], object='ligand')
    ([min_x, min_y, min_z], [max_x, max_y, max_z]) = cmd.get_extent('ligand')

    min_x = min_x - float(extending)
    min_y = min_y - float(extending)
    min_z = min_z - float(extending)
    max_x = max_x + float(extending)
    max_y = max_y + float(extending)
    max_z = max_z + float(extending)

    size_x = max_x - min_x
    size_y = max_y - min_y
    size_z = max_z - min_z
    center_x = (max_x + min_x) / 2
    center_y = (max_y + min_y) / 2
    center_z = (max_z + min_z) / 2

    cmd.delete('all')

    centroid = np.array((center_x, center_y, center_z))
    box_dims = np.array((size_x, size_y, size_z))

    if verbose:
        return {'center_x': center_x, 'center_y': center_y, 'center_z': center_z}, {'size_x': size_x, 'size_y': size_y,
                                                                                    'size_z': size_z}
    else:
        return centroid, box_dims


def fix_protein(filename='', addHs_pH=7.4, output='', try_renumberResidues=True, removeHeterogens=True):
    fix = PDBFixer(filename=filename)
    fix.findMissingResidues()
    fix.findNonstandardResidues()
    fix.replaceNonstandardResidues()
    if removeHeterogens == True:
        fix.removeHeterogens(True)
    fix.findMissingAtoms()
    fix.addMissingAtoms()
    fix.addMissingHydrogens(addHs_pH)
    PDBFile.writeFile(fix.topology, fix.positions, open(output, 'w'))

    if try_renumberResidues == True:
        try:
            original = mda.Universe(filename)
            from_fix = mda.Universe(output)

            resNum = [res.resid for res in original.residues]
            for idx, res in enumerate(from_fix.residues):
                res.resid = resNum[idx]

            save = PDB.PDBWriter(filename=output)
            save.write(from_fix)
            save.close()
        except Exception:
            print('Not possible to renumber residues, check excepton for extra details')


if __name__ == '__main__':
    protein = '/home/anton/PycharmProjects/cadd/data/prots/cox.pdb'
    # centroid, box_dims = get_box(protein, extending=5.0)
    # print(centroid, box_dims, sep='\n')
    fix_protein(filename=protein, addHs_pH=7.4, output='cox_fix.pdb', try_renumberResidues=True, removeHeterogens=False)
