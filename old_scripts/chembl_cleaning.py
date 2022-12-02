import pandas as pd
from rdkit import Chem
from rdkit.Chem.Fragments import fr_Al_COO, fr_Ar_COO


# start: all chembl cox-2 inhibitors with IC50, and already in 3 or 4 phase
df = pd.read_csv('cox_chembl.csv', header=0, error_bad_lines=False, sep=';',
                 usecols=['Molecule ChEMBL ID','Molecule Name',  'Smiles', 'Molecule Max Phase',
                          'Standard Type', 'Standard Relation', 'Standard Value',
       'Standard Units'])
print(df.head())
print('Number of compound in original dataset:', len(df))
# check for nan values
print(df.isna().sum())
# drop compounds with no Standard Value
df.dropna(inplace=True, subset=['Standard Value'])
# print('Number of compound after dropping compounds with no Standard Value:', len(df))
# check for nan values
print(df.isna().sum())
# Sort compounds by Standard Value
df.sort_values(by=['Standard Value'], inplace=True)
# drop duplicates, keep the first one
df.drop_duplicates(subset=['Molecule ChEMBL ID'], inplace=True)
print('Number of compound after dropping duplicates of same Chembl ID:', len(df))
df.reset_index(drop=True, inplace=True)
def convert_smiles_to_mol(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol
    except:
        return None
# convert smiles to mol
df['Mol'] = df['Smiles'].apply(convert_smiles_to_mol)
# drop compounds with no mol
df.dropna(inplace=True, subset=['Mol'])
print('Number of compound after dropping compounds with no mol:', len(df))
df['fr_Al_COO'] = df['Mol'].apply(fr_Al_COO)
df['fr_Ar_COO'] = df['Mol'].apply(fr_Ar_COO)
df = df[df['fr_Al_COO'] + df['fr_Ar_COO'] == 1]
print('Number of compound, which are acids:', len(df))

from commercial_availability.checker import Molport
df['is_commercially_available'] = df['Smiles'].apply(lambda x: Molport().find_compound(smiles=x).is_commercial)
print('Number of compounds, which are commercially available:', len(df[df['is_commercially_available'] == True]))
df.to_csv('cox_chembl_cleaned_v3.csv', index=False)

