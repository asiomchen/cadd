import time
import random
from pathlib import Path
import matplotlib.pyplot as plt
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
from scipy.stats import mode

class Automaton:
    def __init__(self):
        self.data = None
        self.smi = None
        self.compounds = []
        self.fps = []
        self.similarity_matrix = None
        self.clusters = None
        self.triangular_matrix = None

    def load_chembl(self, path, sep=';', sanitize=True) -> None:
        data = pd.read_csv(path, sep=';')
        if sanitize:
            self.data = self._clean_df(data)
        else:
            self.data = data
        self.smi = self.data['Smiles'].tolist()
        self.compounds, self.fps = self._get_fingerprint(self.smi)

    def load_smi(self, path) -> None:
        with open(path, 'r') as f:
            self.smi = f.readlines()

        self.compounds, self.fps = self._get_fps_smi()

    def _get_fingerprint(self, smi: list) -> np.array:
        compounds = []
        for _, chembl_id, smiles in self.data[["Molecule ChEMBL ID", "Smiles"]].itertuples():
            compounds.append((Chem.MolFromSmiles(smiles), chembl_id))

        rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
        fingerprints = [rdkit_gen.GetFingerprint(mol) for mol, idx in compounds]
        return compounds, fingerprints

    def _get_fps_smi(self):
        compounds = [Chem.MolFromSmiles(smi) for smi in self.smi]
        rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
        fingerprints = [rdkit_gen.GetFingerprint(mol) for mol in compounds]
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
    # a.load_chembl('../data/cox-2_chembl.csv')
    # a.butina_cluster()
    # print(a.clusters)
    # print(len(a.clusters))
    a.load_smi('/home/anton/PycharmProjects/cadd/products_long.smi')
    # print(a.smi)
    # print(a.compounds)
    # print(a.fps)
    print(len(set(a.smi)))
    n_clusters = []
    ts = []
    results = {}
    for t in range(0, 35, 1):
        t = t / 100
        ts.append(t)
        a.butina_cluster(t)
        print('Threshold: ', t)
        print(f'Number of clusters: {len(a.clusters)}')
        n_clusters.append(len(a.clusters))
        print(f'Number of compounds: {len(a.compounds)}')
        clusters = a.clusters
        clusters = sorted(clusters, key=len, reverse=True)
        compounds = a.compounds
        fingerprints = a.fps
        print(
            f"Number of clusters: {len(clusters)} from {len(compounds)} molecules at distance cut-off {t:.2f}"
        )
        print("Number of molecules in largest cluster:", len(clusters[0]))
        print("Number of molecules in smallest cluster:", len(clusters[-1]))
        if len(clusters[-1]) == 1:
            print('Singletons: ', sum([1 for cluster in clusters if len(cluster) == 1]))
            print('Doubles: ', sum([1 for cluster in clusters if len(cluster) == 2]))
            print('Triples: ', sum([1 for cluster in clusters if len(cluster) == 3]))
        # calculate percentage of molecules clusters - all clusters expect singletons and doubles
        molecules_clustered = sum([len(cluster) for cluster in clusters if len(cluster) > 2])
        print('Percentage of molecules in clusters: ', (molecules_clustered / len(compounds)))
        print("Average cluster size:", np.mean([len(c) for c in clusters]))
        print("Median cluster size:", np.median([len(c) for c in clusters]))
        print('Mode cluster size:', mode([len(c) for c in clusters]))
        random_idx_1 = np.random.randint(0, len(clusters[0]))
        random_idx_2 = np.random.randint(0, len(clusters[0]))
        random_idx_3 = np.random.randint(0, len(clusters[1]))
        print(
            f"Similarity between two random points in same cluster: {DataStructs.TanimotoSimilarity(fingerprints[clusters[0][random_idx_1]], fingerprints[clusters[0][random_idx_2]]):.2f}"
        )
        print(
            f"Similarity between two random points in different cluster: {DataStructs.TanimotoSimilarity(fingerprints[clusters[0][random_idx_1]], fingerprints[clusters[1][random_idx_3]]):.2f}"
        )
        print('-------------------------------------')
        ## add all printed values to results dict for convertiong to df
        modes = mode([len(c) for c in clusters])
        mode_value = modes[0][0]
        mode_count = modes[1][0]
        results[t] = [len(clusters), len(compounds), len(clusters[0]), len(clusters[-1]),
                      np.mean([len(c) for c in clusters]), np.median([len(c) for c in clusters]), mode_value,
                      mode_count, DataStructs.TanimotoSimilarity(fingerprints[clusters[0][random_idx_1]],
                                                                 fingerprints[clusters[0][random_idx_2]]),
                      DataStructs.TanimotoSimilarity(fingerprints[clusters[0][random_idx_1]],
                                                     fingerprints[clusters[1][random_idx_3]]),
                      (molecules_clustered / len(compounds))]
        print(results[t])
    df = pd.DataFrame.from_dict(results, orient='index', columns=['Number of clusters', 'Number of compounds',
                                                                  'Number of molecules in largest cluster',
                                                                  'Number of molecules in smallest cluster',
                                                                  'Average cluster size', 'Median cluster size',
                                                                  'Mode cluster size', 'Mode cluster count',
                                                                  'Similarity between two random points in same cluster',
                                                                  'Similarity between two random points in different cluster',
                                                                  'Percentage of molecules in clusters'])
    df.to_csv('results.csv')
    plt.plot(ts, n_clusters)
    plt.xlabel('Threshold')
    plt.ylabel('Number of clusters')
    plt.show()
