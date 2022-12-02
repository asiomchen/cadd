import os

from deepchem.utils.rdkit_utils import load_molecule, write_molecule
from dimorphite_dl import DimorphiteDL
from openbabel import pybel


def to_pqbqt(smiles: str, smiles_id: str, silent: bool = False, pH: float = 7.4):
    """
    Convert smiles to pdbqt file for Vina docking
    Args:
        smiles: smiles string
        smiles_id: name of the molecule
        inplace: if True, the converted file will be saved in the same directory as the original file
        silent: if False, exception is thrown when conversion fails. Otherwise, None is returned.

    Returns: None

    """
    try:
        # generate smiles for the protonated ligand
        dimorphite_dl = DimorphiteDL(
            min_ph=pH,
            max_ph=pH,
            max_variants=1,
            label_states=False,
        )
        protonated_smiles = dimorphite_dl.protonate(smiles)[0]
        # creates sdf 3d molecule from smiles using pybel
        sdf_mol = f"{smiles_id}.sdf"
        out = pybel.Outputfile(filename=sdf_mol, format='sdf', overwrite=True)
        mol = pybel.readstring(string=protonated_smiles, format='smiles')
        mol.title = smiles_id
        mol.make3D('mmff94s')
        mol.localopt(forcefield='mmff94s', steps=500)
        out.write(mol)
        out.close()
        # converts the sdf file to pdbqt file using deepchem because opebabel does not support preserve
        # hydrogens on carbon atoms
        output_filename = sdf_mol.replace(".sdf", ".pdbqt")
        ligand_mol = load_molecule(sdf_mol, calc_charges=True, add_hydrogens=False)
        write_molecule(ligand_mol[1], output_filename)
        os.remove(sdf_mol)
    except Exception as e:
        if silent:
            return None
        else:
            raise ValueError(e)
    return output_filename


if __name__ == "__main__":
    os.chdir('/home/anton/in_dev/Docking_tools/master')
    smiles = 'N[C@@H](Cc1ccc(O)cc1)C(O)=O'
    # to_pqbqt(smiles, "is_charged_tyr", pH=7.4)
    dimorphite_dl = DimorphiteDL(
        min_ph=7.4,
        max_ph=7.4,
        max_variants=3,
        label_states=False,
    )
    protonated_smiles = dimorphite_dl.protonate(smiles)
    print(protonated_smiles)