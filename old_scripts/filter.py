import os

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools


def conjugate_ligand(ligand_smiles):
    ligand_rdkit = Chem.MolFromSmiles(ligand_smiles)
    base_ps = 'NCCCCCCNC(=O)c1ccc(cc1)-c1c2ccc(n2)c(-c2c(F)c(F)c(F)c(F)c2F)c2ccc([nH]2)c(-c2c(F)c(F)c(F)c(F)c2F)c2ccc(n2)c(-c2c(F)c(F)c(F)c(F)c2F)c2ccc1[nH]2'
    base_ps_rdkit = Chem.MolFromSmiles(base_ps)
    rxn_smarts = AllChem.ReactionFromSmarts('[CX4:0][NH2:1].[OH:2][C:3]=[0:4]>>[CX4:0][NH1:1][C:3]=[0:4]')
    try:
        product = rxn_smarts.RunReactants((base_ps_rdkit, ligand_rdkit))[0][0]
    except:
        return None

    return Chem.MolToSmiles(product, isomericSmiles=True)

if __name__ == '__main__':
    os.chdir('/home/anton/in_dev/Docking_tools/master')
    df = pd.read_csv('./data/cox-2_chembl_cleaned_acids_filtered.csv')
    df['conjugate'] = df['Smiles'].apply(conjugate_ligand)
    df.drop_duplicates(subset=['Smiles'], inplace=True)
    df.dropna(inplace=True)
    df['Standard Value'] = df['Standard Value'].astype(float)
    df.sort_values(by=['Standard Value'], inplace=True)
    df.to_csv('ligands_conjugated.csv')
