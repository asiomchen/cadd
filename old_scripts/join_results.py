import pandas as pd


# gnina = pd.read_csv('/home/anton/in_dev/Docking_tools/master/gnina_output_conjugates.csv')
# smina = pd.read_csv('/home/anton/in_dev/Docking_tools/master/smina_conjugates.csv')
# # add _smina suffix to smina columns
# smina.columns = ['ID'] + [f'{col}_smina' for col in smina.columns[1:]]
# all = pd.merge(gnina, smina, on='ID', how='outer')
# named = pd.read_csv('/home/anton/in_dev/Docking_tools/master/ligands/cox_chembl_cleaned_v3.csv')
# all = pd.merge(all, named, left_on='ID', right_on='Molecule ChEMBL ID', how='left')
# all.sort_values(by=['Affinity_1_smina'], inplace=True, ascending=True)
# print(all)
# all.to_csv('/home/anton/in_dev/Docking_tools/master/joined_results_smina_named.csv', index=False)

ligands_old = pd.read_csv('/master/ligands/ligands_conjugated.csv')
print(ligands_old)
ligands_old.drop(columns=['Unnamed: 0'], inplace=True)
from commercial_availability.checker import Molport
ligands_old['is_available'] = ligands_old['Smiles'].apply(lambda x: Molport().find_compound(smiles=x).is_commercial)
print(ligands_old['is_available'].value_counts())
ligands_old = ligands_old[ligands_old['is_available']]
ligands_old.to_csv('/home/anton/in_dev/Docking_tools/master/ligands_conjugated_available.csv', index=False)
