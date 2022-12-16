import time
import random
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import Draw
from rdkit.Chem import rdFingerprintGenerator
from molecular_docking.smina_docking import docking_smina_single
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn_extra.cluster import KMedoids


class Automaton:
    def __init__(self):
        self.data = None
        self.smi = None
        self.compounds = []
        self.fps = []
        self.similarity_matrix = None
        self.clusters = None
        self.triangular_matrix = None

    def load_chembl(self, path, sep=';') -> None:
        data = pd.read_csv(path, sep=';')
        self.data = self._clean_df(data)
        self.smi = self.data['Smiles'].tolist()
        self.compounds, self.fps = self._get_fingerprint(self.smi)

    def _get_fingerprint(self, smi: list) -> np.array:
        compounds = []
        for _, chembl_id, smiles in self.data[["Molecule ChEMBL ID", "Smiles"]].itertuples():
            compounds.append((Chem.MolFromSmiles(smiles), chembl_id))

        rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
        fingerprints = [rdkit_gen.GetFingerprint(mol) for mol, idx in compounds]
        return compounds, fingerprints

    def _clean_df(self, df: pd.DataFrame) -> pd.DataFrame:
        df = df[['Molecule ChEMBL ID', 'Smiles', 'Standard Value', 'Standard Units']]
        df['Standard Units'] = df['Standard Units'].apply(lambda x: x if x == 'nM' else np.nan)
        df = df.dropna()
        sorted_df = df.sort_values(by='Standard Value')
        sorted_unique = sorted_df.drop_duplicates(keep='first')
        sorted_unique = sorted_unique.drop_duplicates(subset='Molecule ChEMBL ID', keep='first')
        return sorted_unique

    def cluster(self, n=10) -> None:
        """Cluster fingerprints with KMedoids"""
        self.similarity_matrix = tanimoto_similarity_matrix(self.fps)
        self.clusters = KMedoids(n_clusters=n, metric='precomputed').fit(self.similarity_matrix)
        self.data['cluster'] = self.clusters.labels_

    def butina_cluster(self, threshold=0.5) -> None:
        """Cluster fingerprints with Butina"""
        self.triangular_matrix = tanimoto_distance_matrix(self.fps)
        self.clusters = Butina.ClusterData(self.triangular_matrix, len(self.compounds), threshold, isDistData=True)
        # self.data['cluster'] = self.clusters.labels_


def tanimoto_similarity_matrix(fps):
    """
    Compute the tanimoto similarity matrix of a list of molecules
    Args:
        mols: list of molecules

    Returns: tanimoto similarity matrix

    """
    n = len(fps)
    sim_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            sim_mat[i, j] = sim
            sim_mat[j, i] = sim
    return sim_mat


def tanimoto_distance_matrix(fp_list):
    """Calculate distance matrix for fingerprint list"""
    dissimilarity_matrix = []
    # Notice how we are deliberately skipping the first and last items in the list
    # because we don't need to compare them against themselves
    for i in range(1, len(fp_list)):
        # Compare the current fingerprint against all the previous ones in the list
        similarities = DataStructs.BulkTanimotoSimilarity(fp_list[i], fp_list[:i])
        # Since we need a distance matrix, calculate 1-x for every element in similarity matrix
        dissimilarity_matrix.extend([1 - x for x in similarities])
    return dissimilarity_matrix


if __name__ == '__main__':
    a = Automaton()
    a.load_chembl('../data/cox-2_chembl.csv')
    a.cluster()
    print(a.clusters)
