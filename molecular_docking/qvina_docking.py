from protein_preparation.protein_preparation import *
import subprocess
import os, sys
from utils.convert_alt import to_sdf, to_pdbqt
import pandas as pd


class QvinaWrapper:
    def __init__(self,
                 qvina_path: str = '/home/anton/PycharmProjects/cadd/notebooks/qvina-w',
                 protein_file=None,
                 ligand=None,
                 ligand_type='pdbqt',
                 box_def=None,
                 num_modes=5,
                 exhaustiveness=8,
                 seed=42,
                 out_dir='qvina_output'):
        self.qvina_path = qvina_path
        self.protein_file = protein_file
        self.ligand_file = ligand
        self.out_dir = out_dir
        self.num_modes = num_modes
        self.exhaustiveness = exhaustiveness
        self.seed = seed
        self.ligand_type = ligand_type
        self.box_def = box_def

    def set_protein(self, protein):
        self.protein_file = protein

    def set_ligand(self, ligand):
        self.ligand_file = ligand

    def set_box(self, box_def):
        self.box_def = box_def

    def set_num_modes(self, num_modes):
        self.num_modes = num_modes

    def set_exhaustiveness(self, exhaustiveness):
        self.exhaustiveness = exhaustiveness

    def set_seed(self, seed):
        self.seed = seed

    def set_ligaand_type(self, ligand_type):
        self.ligand_type = ligand_type

    def set_ligand_type(self, ligand_type):
        self.ligand_type = ligand_type

    def set_out_dir(self, out_dir):
        self.out_dir = out_dir

    def _get_box_dims(self):
        if self.box_def is None:
            centroid, size = get_box(self.protein_file)
        else:
            centroid, size = self.box_def
        return centroid, size

    def _pdbqt_from_smiles(self, smiles, smiles_id='ligand'):
        self.ligand_file = to_pdbqt(smiles, smiles_id)

    def docking(self, mol_name='docked', blind=True):
        if blind:
            centroid, size = self._get_box_dims()
            x, y, z = centroid
            x_size, y_size, z_size = size

        self.out_dir = os.path.join(os.getcwd(), self.out_dir)
        if self.ligand_type == 'pdbqt':
            self.ligand_file = self.ligand_file
        elif self.ligand_type == 'smiles':
            self.ligand_file = to_pdbqt(self.ligand_file, mol_name)
        mol_name = mol_name + '.pdbqt'
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        cmd = f'{self.qvina_path} --receptor {self.protein_file} --ligand {self.ligand_file} --num_modes {self.num_modes} --exhaustiveness {self.exhaustiveness} --seed {self.seed} --out {self.out_dir + "/" + mol_name}'
        cmd += f' --center_x {x} --center_y {y} --center_z {z} --size_x {x_size} --size_y {y_size} --size_z {z_size}'
        print(cmd)
        print(self.out_dir)
        subprocess.call(cmd, shell=True)


class QvinaAnalyzer:
    def __init__(self, qvina_output_dir):
        self.qvina_output_dir = qvina_output_dir

    def __call__(self, *args, **kwargs):
        return self.get_scores(*args, **kwargs)

    def get_scores(self, name='docked', as_df=False):
        scores = []
        name += '.pdbqt'
        for file in os.listdir(self.qvina_output_dir):
            if name in file:
                with open(os.path.join(self.qvina_output_dir, file)) as f:
                    for line in f:
                        if 'REMARK VINA RESULT:' in line:
                            scores.append(line.split()[3])
        if as_df:
            scores = pd.DataFrame(scores, columns=['score'])
        return scores


if __name__ == '__main__':
    # protein = '/home/anton/PycharmProjects/cadd/data/prots/cox.pdbqt'
    #
    # ligand = '/home/anton/PycharmProjects/cadd/data/prots/ind_conj.pdbqt'
    # ligand_smi = 'COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1'
    # qvina = QvinaWrapper()
    # qvina.set_protein(protein)
    # qvina.set_ligand(ligand_smi)
    # qvina.set_ligand_type('smiles')
    # qvina.set_num_modes(10)
    # qvina.set_exhaustiveness(64)
    # qvina.set_seed(42)
    # qvina.set_out_dir('qvina_output')
    # qvina.docking('docked_64_smi')
    analyzer = QvinaAnalyzer('qvina_output')
    print(analyzer('docked_64_smi', as_df=False))
