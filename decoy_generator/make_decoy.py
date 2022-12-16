import pandas as pd
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Chem import DataStructs, rdMolDescriptors, Crippen

from fixes.convert_alt import to_rdkit_mol


# based on approach from https://academic.oup.com/bioinformatics/article/28/12/1661/270210

def pandas_add_descriptors(df: pd.DataFrame, fingerprint='Morgan') -> pd.DataFrame:
    '''
    Add descriptors to a dataframe - Molecular Weight, Number of Rotational Bonds,
    Number of HBDs, LogP, Fingerprint(MACCS, Morgan)
    :param df:
    :return: modified dataframe
    '''
    # Add the molecular weight
    df['molecular_weight'] = df['mol'].apply(lambda x: rdMolDescriptors.CalcExactMolWt(x))
    # Add the number of rotational bonds
    df['rotational_bonds'] = df['mol'].apply(lambda x: rdMolDescriptors.CalcNumRotatableBonds(x))
    # Add the number of HBDs
    df['hbd'] = df['mol'].apply(lambda x: rdMolDescriptors.CalcNumHBD(x))
    # Add the log P-value
    df['log_p'] = df['mol'].apply(lambda x: Crippen.MolLogP(x))
    if fingerprint == 'Morgan':
        # Add the Morgan fingerprint
        df['fingerprint'] = df['mol'].apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 3, 1024))
    if fingerprint == 'MACCS':
        # Add the MACCS keys
        df['fingerprint'] = df['mol'].apply(lambda x: MACCSkeys.GenMACCSKeys(x))
    return df

def compare_mols(mol, decoy_mol, similarity_threshold=0.75)-> bool:
    '''
    Compare the molecular weight, number of rotational bonds,
    number of HBDs, and logP value of two molecules. Tanimoto similarity is used to compare the fingerprints.

    :param mol: rdkit molecule
    :param decoy_mol: rdkit molecule
    :param similarity_threshold: threshold for Tanimoto similarity
    :return: True if decoy_mol is decoy of mol
    '''
    # Get the molecular weight
    molecular_weight = rdMolDescriptors.CalcExactMolWt(mol)
    decoy_molecular_weight = rdMolDescriptors.CalcExactMolWt(decoy_mol)
    if abs(molecular_weight - decoy_molecular_weight) > 25:
        return False

    # Get the number of rotational bonds
    rotational_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    decoy_rotational_bonds = rdMolDescriptors.CalcNumRotatableBonds(decoy_mol)
    if abs(rotational_bonds - decoy_rotational_bonds) > 1:
        return False

    # Get the number of HBDs
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    decoy_hbd = rdMolDescriptors.CalcNumHBD(decoy_mol)
    if abs(hbd - decoy_hbd) > 2:
        return False

    # Get the log P-value
    log_p = Crippen.MolLogP(mol)
    decoy_log_p = Crippen.MolLogP(decoy_mol)
    if abs(log_p - decoy_log_p) > 1.0:
        return False
    # Get the Tanimoto similarity of the fingerprints
    fingerprint_mol = MACCSkeys.GenMACCSKeys(mol)
    fingerprint_decoy_mol = MACCSkeys.GenMACCSKeys(decoy_mol)
    if DataStructs.TanimotoSimilarity(fingerprint_mol, fingerprint_decoy_mol) > similarity_threshold:
        return False
    return True


def check_new_decoy(mol, decoy_mols, similarity_threshold=0.9)-> bool:
    '''
    Check if a new decoy is different from any of the previously selected decoys
    :param mol: smiles of the potential decoy
    :param decoy_mols: list of smiles of the previously selected decoys
    :return: True if the new decoy is different from any of the previously selected decoys
    '''
    mol = to_rdkit_mol(mol)
    mol_fingerprint = MACCSkeys.GenMACCSKeys(mol)
    for existing_decoy in decoy_mols:
        existing_decoy = to_rdkit_mol(existing_decoy)
        decoy_fingerprint = MACCSkeys.GenMACCSKeys(existing_decoy)
        if DataStructs.TanimotoSimilarity(mol_fingerprint, decoy_fingerprint) > similarity_threshold:
            return False
    return True


def run_search(csv_file='chembl_dataset.smi.txt', mol=None, mol_rdkit=None, n=5) -> list:
    '''
    Run the decoy search algorithm
    :param csv_file:  path to the csv file containing the dataset
    :param mol:  smiles of the target molecule
    :param mol_rdkit:  rdkit molecule of the target molecule
    :param n: Number of decoys to be generated
    :return:
    '''
    with pd.read_csv(csv_file, chunksize=50000, delimiter=' ') as reader:
        end = False
        out = []
        for chunk in reader:
            chunk = chunk.sample(frac=1)
            chunk.columns = ['smiles', 'id']
            chunk['mol'] = chunk['smiles'].apply(lambda x: to_rdkit_mol(x))
            # chunk = pandas_add_descriptors(chunk, fingerprint='MACCS')
            if end:
                break
            for smiles, decoy_mol in zip(chunk['smiles'], chunk['mol']):
                if compare_mols(mol_rdkit, decoy_mol):
                    if len(out) != 0:
                        if check_new_decoy(mol, out):
                            out.append(smiles)

                    else:
                        out.append(smiles)
                    if len(out) == n:
                        return out
                    else:
                        print(f'Found {len(out) + 1} decoys')



def find_decoys(mol, decoy_source='chembl', mol_type='smiles', n=5) -> list:
    results = []
    print('Any decoy source can be used')
    sources = ['chembl']
    if mol_type == 'smiles':
        mol_rdkit = to_rdkit_mol(mol)
    if decoy_source not in sources:
        print('Decoy source not supported')
        return

    if decoy_source == 'chembl':
        # load source in chunks of 100 and descriptors for each chunk and save to csv
        results = run_search(csv_file='chembl_dataset.smi.txt', mol=mol, mol_rdkit=mol_rdkit, n=n)
    return results



if __name__ == '__main__':
    print(find_decoys('Cc1c(cccc1c2ccccc2)COc3ccc(c(n3)OC)CNCCNC(=O)C', decoy_source='chembl', mol_type='smiles', n=1))





