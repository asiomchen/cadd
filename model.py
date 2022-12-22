import numpy as np
from rdkit.DataStructs import BulkTanimotoSimilarity

from create_conjugates.reaction import Reactor
from rdkit import Chem
from itertools import product, permutations
from rdkit.Chem import DataStructs
from rdkit import Chem


class Model:
    def __init__(self, ps_path, linker_path, ligand_path):
        self.ps_path = ps_path
        self.linker_path = linker_path
        self.ligand_path = ligand_path
        self.pss = open(ps_path, 'r').readlines()
        self.linkers = open(linker_path, 'r').readlines()
        self.ligands = open(ligand_path, 'r').readlines()
        self.reactor = Reactor()
        self.linkers_space = [Chem.MolFromSmiles(x.strip()) for x in self.linkers]
        self.pss_space = [Chem.MolFromSmiles(x.strip()) for x in self.pss]
        self.ligands_space = [Chem.MolFromSmiles(x.strip()) for x in self.ligands]
        self.compounds_space = self.get_compounds_space()
        # self.linkers_similarity = self.get_similarity('linkers')
        # self.ps_similarity = self.get_similarity('ps')
        # self.ligands_similarity = self.get_similarity('ligands')

    def get_compounds_space(self):
        '''Returns an array of all possible compounds'''
        reagents_space = list(product(self.pss_space, self.linkers_space, self.ligands_space))
        print(reagents_space)
        print(f'Number of ps: {len(self.pss_space)}')
        print(f'Number of linkers: {len(self.linkers_space)}')
        print(f'Number of ligands: {len(self.ligands_space)}')
        print(self.ligands)
        print(f'Number of possible compounds: {len(list(reagents_space))}')
        reactor = Reactor()
        compounds_space = [reactor.conjugate(*x) for x in reagents_space]
        print(f'Number of compounds generated: {len(compounds_space)}')
        return compounds_space

    def get_similarity(self, name):
        '''Returns an array of all possible compounds'''
        if name == 'linkers':
            return self.get_similarity_from_list(self.linkers)
        elif name == 'ps':
            return self.get_similarity_from_list(self.pss)
        elif name == 'ligands':
            return self.get_similarity_from_list(self.ligands)
        else:
            raise ValueError('Wrong name')

    def get_similarity_from_list(self, list_):
        '''Returns an array of all possible compounds'''
        similarity = []
        for i, j in permutations(list_, 2):
            i = Chem.MolFromSmiles(i.strip())
            j = Chem.MolFromSmiles(j.strip())
            similarity.append(self.compare_mol(i, j))
        return similarity

    def compare_mol(self, mol1, mol2):
        fp1 = Chem.RDKFingerprint(mol1)
        fp2 = Chem.RDKFingerprint(mol2)
        return DataStructs.FingerprintSimilarity(fp1, fp2)

    def bulk_predict_similarity(self, cmpd, compounds_space):
        compounds_space_fps = [Chem.RDKFingerprint(x) for x in compounds_space]
        cmpd_fp = Chem.RDKFingerprint(cmpd)
        return BulkTanimotoSimilarity(cmpd_fp, compounds_space_fps)

    def save_products(self, name='products.smi'):
        with open(name, 'w') as f:
            for i in self.compounds_space:
                f.write(Chem.MolToSmiles(i) + '\n')


test_conj = 'COc1ccc2c(c1)c(CC(=O)OCCOCCOCCOCCOC(=O)c1ccc(-c3c4nc(c(-c5cc(S(=O)(=O)[O-])ccc5F)c5ccc([nH]5)c(-c5cc(S(=O)(=O)[O-])ccc5F)c5nc(c(-c6cc(S(=O)(=O)[O-])ccc6F)c6ccc3[nH]6)C=C5)C=C4)cc1)c(C)n2C(=O)c1ccc(Cl)cc1'
test_link = 'OCCOCCO'
agent = Model('./data/parts/ps.smi', './data/parts/linker.smi', './data/parts/ligands.smi')
agent.save_products('products_long.smi')
