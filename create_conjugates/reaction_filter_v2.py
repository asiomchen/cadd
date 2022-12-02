import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import rdDepictor

from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from rdkit.Chem import PandasTools

rdDepictor.SetPreferCoordGen(True)


def create_conjugate(porphyrin: Chem.Mol, ligand: Chem.Mol, target_base='amine', target_ligand='acid',
                     base_ligand_numbers=1) -> Chem.Mol or None:
    """
    Creates a conjugate of a given base porphyrin for a given ligand.
    Args:
        porphyrin: <rdkit.Chem.rdchem.Mol object> - target porphyrin to make conjugate from as Mol object
        ligand: <rdkit.Chem.rdchem.Mol object> - ligand to conjugate with as Mol object
        target_base: str target base group for conjugation - linking group type of the target porphyrin
        target_ligand: str target ligand group for conjugation

    Returns: <rdkit.Chem.rdchem.Mol object> of the conjugate

    """
    # base_ligand_numbers = fr_Al_COO(porphyrin) + fr_Ar_COO(porphyrin)
    # reaction def for basic replacement od acidic proton to At atom
    # to avoid reaction with other -OH groups, than in base porphyrin
    rxn_smarts = '[OH:1]>>[O:1][At]'
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    sub = porphyrin

    # reaction definition for amid formation
    if target_base == 'acid' and target_ligand == 'primary amine':
        rxn_smarts = '[*:0][O:1][At].[*:3][NH2:4]>>[*:0][NH:4][*:3]'
    if target_base == 'amine' and target_ligand == 'acid':
        rxn_smarts = '[*:0][NH2:1].[CX3](=O)[OX2H1]>>[*:0][NH:1][CX3](=O)'
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    # runs reaction in all possible places as defined by base_ligand_numbers

    product = rxn.RunReactants((sub, ligand))[0][0]
    Chem.SanitizeMol(product)

    return product


def conjugate_ligand(ligand_smiles):
    base_ps = 'NCCCCCCNC(=O)c1ccc(cc1)-c1c2ccc(n2)c(-c2c(F)c(F)c(F)c(F)c2F)c2ccc([nH]2)c(-c2c(F)c(F)c(F)c(F)c2F)c2ccc(n2)c(-c2c(F)c(F)c(F)c(F)c2F)c2ccc1[nH]2'
    rxn_smarts = AllChem.ReactionFromSmarts('[CX4:0][NH2:1].[OH:2][C:3]=[0:4]>>[CX4:0][NH1:1][C:3]=[0:4]')
    # rxn_smarts = AllChem.ReactionFromSmarts('[NH2:1]>>[NH1:1]')
    product = rxn_smarts.RunReactants((etn, acetic))[0][0]


if __name__ == '__main__':
    etn = 'NCCCCCCNC(=O)c1ccc(cc1)-c1c2ccc(n2)c(-c2c(F)c(F)c(F)c(F)c2F)c2ccc([nH]2)c(-c2c(F)c(F)c(F)c(F)c2F)c2ccc(n2)c(-c2c(F)c(F)c(F)c(F)c2F)c2ccc1[nH]2'
    acetic = 'CCCCCCC(Sc1nc(Cl)cc(Oc2ccc3ncccc3c2)n1)C(=O)O'

    etn = Chem.MolFromSmiles(etn)
    acetic = Chem.MolFromSmiles(acetic)
    # conj = create_conjugate(pfp, rcx)
    rxn_smarts = AllChem.ReactionFromSmarts('[CX4:0][NH2:1].[OH:2][C:3]=[0:4]>>[CX4:0][NH1:1][C:3]=[0:4]')
    # rxn_smarts = AllChem.ReactionFromSmarts('[NH2:1]>>[NH1:1]')
    product = rxn_smarts.RunReactants((etn, acetic))[0][0]
    df = pd.read_csv('/master/ligands/cox_chembl_cleaned_v3.csv')
    PandasTools.AddMoleculeColumnToFrame(df, smilesCol='Smiles', molCol='Mol')
    df['Conj'] = df['Mol'].apply(lambda x: rxn_smarts.RunReactants((etn, x))[0][0])
    df['Conj_smiles'] = df['Conj'].apply(lambda x: Chem.MolToSmiles(x))
    df.drop(columns=['Mol', 'Conj'], inplace=True)
    df.to_csv('/home/anton/in_dev/Docking_tools/master/cox_chembl_cleaned_v3_conj.csv')
