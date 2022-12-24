from rdkit.Chem.MolStandardize import rdMolStandardize

from create_conjugates.reaction import Reactor
import os, sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, SaltRemover
from utils.convert_alt import to_pdbqt

df = pd.read_csv('pocket_random_subset.csv', sep=',')

ps_smi = r'OC(C1=CC=C(C=C1)/C2=C3C=CC(/C(C4=C(F)C(F)=C(F)C(F)=C4F)=C5N/C(C=C\5)=C(C6=N/C(C=C6)=C(C7=CC=C2N7)/C8=C(F)C(F)=C(F)C(F)=C8F)/C9=C(F)C(F)=C(F)C(F)=C9F)=N\3)=O'
ps = Chem.MolFromSmiles(ps_smi)
link_smi = r'NCCCCCCN'
link = Chem.MolFromSmiles(link_smi)

reactor = Reactor()
conjugates = []
errors = []
for smiles in df['Smiles']:
    smiles = smiles.split('.')[0].replace('[O-]', 'O')
    mol = Chem.MolFromSmiles(smiles)
    try:
        mol = reactor.conjugate(ps, link, mol)
        Chem.SanitizeMol(mol)
        print(Chem.MolToSmiles(mol))
        conjugates.append(Chem.MolToSmiles(mol) + '\n')
    except:
        conjugates.append('None\n')
        errors.append(smiles)
print('Errors: ')
for smiles in errors:
    print(smiles)

with open('pocket_conjugates.smi', 'w') as f:
    f.writelines(conjugates)

ids = df['ID'].tolist()
ps_id = 1
link_id = 1
os.mkdir(f'../notebooks/pocket_set')
os.chdir(f'../notebooks/pocket_set')
for ids, conjugate in zip(ids, conjugates):
    name = f'{ps_id}_{link_id}_{ids}'
    print(conjugate)
    print(name)
    if conjugate != 'None':
        to_pdbqt(conjugate, name)
