from vina import Vina
import logging
import os
import subprocess
import warnings
from typing import List

import pandas as pd
from deepchem.utils.rdkit_utils import load_molecule, write_molecule
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from vina import Vina

logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore")


def optimize_conformation(mol):
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def _exec_subprocess(command: List[str], timeout: int = None) -> List[str]:
    cmd = ' '.join([str(entry) for entry in command])

    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, timeout=timeout)
        out, err, return_code = str(result.stdout, 'utf-8').split('\n'), str(result.stderr, 'utf-8'), result.returncode

        if return_code != 0:
            logger.error('Docking failed with command "' + cmd + '", stderr: ' + err)
            raise ValueError('Docking failed')

        return out
    except subprocess.TimeoutExpired:
        logger.error('Docking failed with command ' + cmd)
        raise ValueError('Docking timeout')


def docking_vina_alt(protein_file: str,
                     smiles_ligand: str,
                     smiles_id: str,
                     num: int,
                     box_def: tuple = (
                             11.338500022888184, -19.122999668121338, 181.1580047607422, 29.503000259399414,
                             18.498000144958496,
                             17.98199462890625),
                     out_dir: str = "vina_output",
                     timeout: int = 600):
    # creation of the output directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # creation of the output file is no longer works, because the output .pdbqt file
    # do not contain all explicit hydrogens only hydrogen on N,O,S atoms is preserved
    # TODO: make function smiles_to_pdbqt() to create new valid pdbqt input file
    # ligand preparation
    mol = Chem.MolFromSmiles(smiles_ligand)
    mol = optimize_conformation(mol)
    Chem.MolToMolFile(mol, 'molecule.mol')

    cmd = 'obabel --imol molecule.mol --osdf ligand.sdf'
    subprocess.run(cmd, shell=True)
    os.remove('molecule.mol')

    ligand_mol = load_molecule("ligand.sdf",
                               calc_charges=True,
                               add_hydrogens=True)
    write_molecule(ligand_mol[0], f"{smiles_id}.pdbqt")

    out_name = f"{smiles_id}.pdbqt"
    out = pybel.Outputfile(filename=out_name, format='pdbqt', overwrite=True)
    mol = pybel.readstring(string=smiles_ligand, format='smiles')
    mol.title = smiles_id
    mol.make3D('mmff94s')
    mol.localopt(forcefield='mmff94s', steps=1000)
    out.write(mol)
    out.close()
    # ligand preparation end

    # vina docking

    v = Vina(sf_name='vina')
    # seting up the receptor
    v.set_receptor(protein_file)
    # setting up the ligand
    v.set_ligand_from_file(out_name)

    v.compute_vina_maps(center=[box_def[0], box_def[1], box_def[2]],
                        box_size=[box_def[3], box_def[4], box_def[5]])
    # Dock the ligand
    # TODO: make exhaustiveness and n_poses to be configurable
    v.dock(exhaustiveness=8, n_poses=1)

    # seting up the output filename
    docked_out_name = f"{num}-{smiles_id}_docked.pdbqt"

    # get current directory
    cur_dir = os.getcwd()

    # seting up the output directory
    os.chdir(out_dir)

    # writing the docked ligand results to the output directory
    v.write_poses(docked_out_name, n_poses=1, overwrite=True)

    # return to the current directory
    os.chdir(cur_dir)

    # removing the ligand input file
    os.remove(out_name)


def drop_error(
        df: pd.DataFrame,
        chembl_ID: str,
        source="random_matrix_results_conjugated_v2_single.csv"
) -> pd.DataFrame:
    """
    Drop the rows with error in the docking process and updating source csv file
    Args:
        df:  DataFrame with the docking results
        chembl_ID: ID of the molecule to be removed
        source: source csv file path

    Returns:
        DataFrame with the docking results without the error rows

    """
    # removing the error row
    df = df.drop(df[df.chembl_ID == chembl_ID].index)
    # updating the source csv file
    df.to_csv(source, index=False)
    return df


if __name__ == "__main__":
    os.chdir('/home/anton/in_dev/Docking_tools/master')
    docking_vina_alt(protein_file="data/A_B_1cx2.pdbqt",
                     smiles_ligand='[C]1=C(N(N=C1C(F)(F)F)C1=[C][C]=C([C]=[C]1)S(=O)(=O)N)C1=[C][C]=C([C]=[C]1)Br',
                     smiles_id='original',
                     num=2)
