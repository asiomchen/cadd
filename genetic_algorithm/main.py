import os, sys

sys.path.insert(0, '..')
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
from utils.convert_alt import to_pdbqt, to_sdf
from utils.fixes import global_seed
from utils.extract_scores import extract_scores
from itertools import product
import numpy as np
from create_conjugates.reaction import Reactor
import random
from rdkit import Chem
import subprocess
import argparse
from multiprocessing import Pool
from other import retrive_chkpt

global_seed()

parser = argparse.ArgumentParser()
# add arguments: population size, number of generations, elite size, parent size
parser.add_argument('--population_size', type=int, default=50, help='Population size', required=False)
parser.add_argument('--generations', type=int, default=20, help='Number of generations', required=False)
parser.add_argument('--elite_size', type=int, default=0.1, help='Elite size', required=False)
parser.add_argument('--parent_size', type=int, default=0.5, help='Parent size', required=False)
parser.add_argument('--mutation_rate', type=float, default=0.333, help='Mutation rate', required=False)
parser.add_argument('--output', type=str, help='Output file', default='output.log')
parser.add_argument('--chkpt', type=str, help='Checkpoint file', default='./docked_mols', required=False)
parser.add_argument('--n_jobs', type=int, help='Number of jobs', default=2, required=False)
# end, parse arguments
args = parser.parse_args()

data_path = '../data/parts/'
ligs = open(data_path + 'ligands_v2.smi').readlines()
links = open(data_path + 'linker.smi').readlines()
pss = open(data_path + 'ps.smi').readlines()
all_compounds = list(product(ligs, links, pss))
print(len(all_compounds))
SCORES_DF = None
if args.chkpt and 'csv' in args.chkpt:
    SCORES_DF = pd.read_csv(args.chkpt, index_col=0)
    print(f'Loaded {len(SCORES_DF)} molecules from {args.chkpt}')
elif args.chkpt and './' in args.chkpt:
    SCORES_DF = retrive_chkpt(args.chkpt)
    print(f'Loaded {len(SCORES_DF)} molecules from folder {args.chkpt}')


def global_sanitize(smiles: str) -> str:
    '''
    Simple function to sanitize smiles, by hydrogen normalization and removing salts
    :param smiles:
    :return:
    '''
    # Remove salts
    smiles = smiles.split('.')[0]
    # Remove hydrogen isotopes
    smiles = smiles.replace('[2H]', '[H]')
    smiles = smiles.replace('[3H]', '[H]')
    return smiles


ligs = list(map(lambda x: global_sanitize(x), ligs))
links = list(map(lambda x: global_sanitize(x), links))
pss = list(map(lambda x: global_sanitize(x), pss))


def MolWithoutIsotopesToSmiles(mol: Chem.Mol) -> str:
    """
    Convert a rdkit molecule to smiles without isotopes
    :param mol:
    :return:
    """
    atom_data = [(atom, atom.GetIsotope()) for atom in mol.GetAtoms()]
    for atom, isotope in atom_data:
        if isotope:
            atom.SetIsotope(0)
    smiles = Chem.MolToSmiles(mol)
    for atom, isotope in atom_data:
        if isotope:
            atom.SetIsotope(isotope)
    return smiles


class Compound:
    '''
    Base class for fragment-base compound descriptor.
    '''

    def __init__(self, ps, link, lig, name='compound', parent_1=None, parent_2=None, is_mutated=False):
        self.ps = self._sanitize(ps)
        self.link = self._sanitize(link)
        self.lig = self._sanitize(lig)
        self._conjugate = None
        if self.conjugate is None:
            print('Error in reaction for compound: ', name)
            print('Scoring will be skipped')
        self.name = name
        self.is_generated = os.path.exists('./generated_mols/' + self.name + '.pdbqt')
        self.is_docked = os.path.exists('./docked_mols/' + self.name + '.pdbqt')
        self.is_mutated = is_mutated
        self.score = 0
        self.parent_1 = parent_1
        self.parent_2 = parent_2

    @property
    def conjugate(self):
        return self._conjugate

    @conjugate.setter
    def conjugate(self, value):
        self._conjugate = value

    @conjugate.getter
    def conjugate(self):
        if self._conjugate is None:
            self._conjugate = self._reaction()
        return self._conjugate

    def __str__(self):
        return f'Compound(name={self.name}, parent_1={self.parent_1.name if self.parent_1 is not None else None}, parent_2={self.parent_2.name if self.parent_2 is not None else None}, is_mutated={self.is_mutated})'

    def __repr__(self):
        return f'Compound({self.name})'

    def __eq__(self, other):
        return self.ps == other.ps and self.link == other.link and self.lig == other.lig and False

    def __hash__(self):
        return hash(self.conjugate)

    def reset(self):
        '''
        Reset generated and docked flags and apply sanitization
        :return:
        '''
        self.is_generated = os.path.exists('./generated_mols/' + self.name + '.pdbqt')
        self.is_docked = os.path.exists('./docked_mols/' + self.name + '.pdbqt')
        self.ps = self._sanitize(self.ps)
        self.link = self._sanitize(self.link)
        self.lig = self._sanitize(self.lig)


    def dock(self, protein='../data/prots/cox2.pdbqt', ex=32, centroid=(42.84, 31.02, 32.31, 34, 75, 43.79, 34.82),
             out_dir='./docked_mols'):
        scores = []
        if self.conjugate is None:
            self.score = 0
            return self.score
        ligand = './generated_mols/' + self.name + '.pdbqt'
        if ligand is None:
            print('Error in generation, smiles: ', self.conjugate)
            return None
        else:
            ligand = './generated_mols/' + self.name + '.pdbqt'

        cmd_template = "./qvina-w --receptor {} --ligand {} --num_modes 3 " \
                       "--exhaustiveness {} --seed 42 --out {} " \
                       "--center_x 42.8405 --center_y 31.0155 --center_z 32.3135 " \
                       "--size_x 34.751 --size_y 43.7889 --size_z 34.821"
        out_name = out_dir + '/' + self.name + '.pdbqt'
        cmd = cmd_template.format(protein, ligand, ex, out_name)
        print(f'Entering docking for {self.name}')
        return_code = subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if return_code != 0:
            print('Error in docking, smiles: ', self.conjugate)
            return 0
        else:
            scores = extract_scores(out_name)
            print(f'Exiting docking for {self.name}, scores: {scores}')
            return np.mean(scores)

    def _sanitize(self, smiles: str) -> str:
        '''
        Simple function to sanitize smiles, by hydrogen normalization and removing salts
        :param smiles:
        :return:
        '''
        smiles = smiles.split('.')[0]
        smiles = smiles.replace('[2H]', '[H]')
        smiles = smiles.replace('[3H]', '[H]')
        return smiles

    def _reaction(self):
        '''
        Perform reaction between ps, link and lig using Reactor class
        Only carboxylic acid-PS, diamine dialcohol linkers and carboxylic acid-ligands reactions are supported
        :return:
        '''
        smiles_compounds = [self.ps, self.link, self.lig]
        # convert smiles to rdkit molecules
        mols = list(map(Chem.MolFromSmiles, smiles_compounds))
        reactor = Reactor()
        # perform reaction
        mol = reactor.conjugate(*mols)
        # convert rdkit molecule to smiles
        if mol is None:
            return mol
        # sanitize Mol
        Chem.SanitizeMol(mol)
        # smiles is returned
        smiles = Chem.MolToSmiles(mol)
        return smiles

    def cross(self, other, mutation_rate=0.333):
        '''
        Crossover of two compounds (conjugates) with probability of mutation
        :param other: Compound object
        :param mutation_rate: float - probability of mutation
        :return: Compound object
        '''
        # probability of each parental chromosome to be chosen
        p = [0.5, 0.5]
        # generate mutation mask
        mask = np.random.choice([0, 1], size=3, p=[1 - mutation_rate, mutation_rate])
        # check if mutation is needed
        is_mutated = any(mask)
        # choose parent for ps and perform mutation if needed
        ps = np.random.choice([self.ps, other.ps], p=p)
        if mask[0] == 1:
            ps = np.random.choice(pss)
        # choose parent for link and perform mutation if needed
        link = np.random.choice([self.link, other.link], p=p)
        if mask[1] == 1:
            link = np.random.choice(links)
        # choose parent for lig and perform mutation if needed
        lig = np.random.choice([self.lig, other.lig], p=p)
        if mask[2] == 1:
            lig = np.random.choice(ligs)
        # name generation based on id of fragments
        name = f'{pss.index(ps)}_{links.index(link)}_{ligs.index(lig)}'
        # create new Compound object
        offspring = Compound(ps, link, lig, name=name, parent_1=self, parent_2=other, is_mutated=is_mutated)
        # if conjugate is not formed, return ERROR Compound, which will have score 0
        if offspring.conjugate is None:
            print('Error in reaction for compound in .cross(): ', name)
            print('Scoring will be skipped')
            offspring.name = 'ERROR_' + name
        # reset is used to ensure that generated and docked flags are set correctly
        # (rather buggy and needs to be fixed with @property for example)
        offspring.reset()
        return offspring

    def mutate(self, mutation_rate=0.1):
        '''
        Mutation of a compound (conjugate) with probability of mutation - legacy function, was used in previous version,
        but currently mutation is performed in .cross() function
        :param mutation_rate:
        :return:
        '''
        new_compound = Compound(self.ps, self.link, self.lig, name=self.name, parent_1=self)
        if random.random() < mutation_rate:
            new_ps = random.choice(pss)
            new_compound.ps = new_ps
        if random.random() < mutation_rate:
            new_link = random.choice(links)
            new_compound.link = new_link
        if random.random() < mutation_rate:
            new_lig = random.choice(ligs)
            new_compound.lig = new_lig
        new_compound.conjugate = new_compound._reaction()
        if new_compound.conjugate is None:
            print('Error in reaction for compound in mutate() method: ', new_compound.name)
            print('Scoring will be skipped')
            new_compound.name = self.name + f'_ERROR_{self.ps}_{self.link}_{self.lig}'
        else:
            mutant_name = f'{pss.index(new_compound.ps)}_{links.index(new_compound.link)}_{ligs.index(new_compound.lig)}'
            new_compound.name = mutant_name
        new_compound.reset()
        return new_compound

    def to_pdbqt(self, path='./generated_mols'):
        '''
        Convert smiles to pdbqt with meeko
        :param path:
        :return:
        '''
        if self.conjugate is None:
            return None
        os.chdir(path)
        # firstly, smiles is converted to sdf with openbabel and protonation state is set to pH 7.4
        sdf_name = to_sdf(self.conjugate, self.name)
        # then, sdf is converted to pdbqt with meeko, as it is preferred way to generate pdbqt files
        # to avoid connectivity and atom type issues
        meeko_cmd = f'mk_prepare_ligand.py -i {sdf_name} -o {self.name}.pdbqt'
        return_code = subprocess.call(meeko_cmd, shell=True)
        # remove sdf file if exists ????
        if os.path.exists(sdf_name):
            os.remove(sdf_name)
        os.chdir('../')
        # return None if conversion was not successful
        if return_code != 0:
            print('Error in meeko, smiles: ', self.conjugate)
            return None
        else:
            # set generated flag to True
            self.is_generated = True
            return self.name + '.pdbqt'


class GeneticDocker:
    def __init__(self, population_size):
        self.dock_calls = 0
        self.ligs = None
        self.links = None
        self.pss = None
        self.population_size = population_size
        self.history = {}

        self.current_iteration = 0
        if not any([self.ligs, self.links, self.pss]) is None:
            self.current_population = None
        else:
            self.current_population = self._init_population(size=self.population_size)
        self.next_population = None
        self.init_genes = None

    def genes_overview(self, population=None):
        '''
        Return proportions of genes in the generation
        '''
        if population is None:
            population = self.current_population
        names_spited = [c.name.split('_') for c in population]
        ps = [c[0] for c in names_spited]
        link = [c[1] for c in names_spited]
        lig = [c[2] for c in names_spited]
        # create dataframe
        df = pd.DataFrame({'ps': ps, 'link': link, 'lig': lig})
        return df

    def _init_population(self, size):
        init_population = []
        for i in range(size):
            ps = random.sample(self.pss, 1)[0]
            link = random.sample(self.links, 1)[0]
            lig = random.sample(self.ligs, 1)[0]
            compound_name = f'{self.pss.index(ps)}_{self.links.index(link)}_{self.ligs.index(lig)}'

            init_population.append(Compound(ps, link, lig, name=compound_name))
        self._generate_population(init_population)
        self.init_genes = self.genes_overview(init_population)
        return init_population

    def _generate_population(self, population):
        comp_count = 0
        to_generate = [comp for comp in population if comp.name + '.pdbqt' not in os.listdir('./generated_mols')]
        comp_count = len(to_generate)
        with Pool(processes=8) as pool:
            pool.map(Compound.to_pdbqt, to_generate)
        print(f'Generated {comp_count} compounds  out of {len(population)}.')

    def set_parts(self, pss, links, ligs):
        self.pss = pss
        self.links = links
        self.ligs = ligs
        self.current_population = self._init_population(size=self.population_size)

    def run_iteration(self):
        global SCORES_DF
        mean_scores = {}
        print('Running iteration: ', self.current_iteration)

        if self.next_population is not None:
            self.current_population = self.next_population
            self.next_population = None
        docking_queue = []
        if SCORES_DF is None:
            print(f'Scores dataframe is empty, running docking for {len(self.current_population)} compounds')
        else:
            pass
            # no_doc = [c for c in self.current_population if c.name not in SCORES_DF.index]
            # print(f'Scores dataframe is not empty, running docking for {len(no_doc)} compounds')
        for compound in self.current_population:
            # print('Docking compound: ', compound.name)
            # print('Compound is docked: ', compound.is_docked)
            if SCORES_DF is None:
                print('SCORES_DF is None')
            if SCORES_DF is not None and compound.name in SCORES_DF.index:
                mean_scores[compound] = SCORES_DF.loc[compound.name, 'score']
            if SCORES_DF is not None and compound.name not in SCORES_DF.index:
                docking_queue.append(compound)
        print(f'Queue length: {len(docking_queue)}')
        print(f'Queue: {docking_queue}')
        with Pool(args.n_jobs) as p:
            # use .dock() method to run docking
            scores = p.map(Compound.dock, docking_queue)
        print(f'Scores length: {len(scores)}')
        new_scores = {compound.name: score for compound, score in zip(docking_queue, scores)}
        print(f'New scores: \n{new_scores}')
        for compound, score in zip(docking_queue, scores):
            mean_scores[compound] = score
            new_scores_df = pd.DataFrame.from_dict(new_scores, orient='index', columns=['score'])
            # add to SCORES_DF
            SCORES_DF = pd.concat([SCORES_DF, new_scores_df])
            # drop index duplicates
            SCORES_DF = SCORES_DF[~SCORES_DF.index.duplicated(keep='first')]

        # apply elitism select only population_size best compounds
        if SCORES_DF is None:
            to_df = {c.name: score for c, score in mean_scores.items()}
            SCORES_DF = pd.DataFrame.from_dict(to_df, orient='index')
            SCORES_DF.columns = ['score']

        self.current_iteration += 1
        # sort by mean score
        sorted_scores = sorted(mean_scores.items(), key=lambda x: x[1])
        print('Applying elitism...')
        self.current_population = [x[0] for x in sorted_scores[:self.population_size]]
        sorted_scores = [x for x in sorted_scores[:self.population_size]]
        print('Current population size: ', len(self.current_population))
        # select 3 mols with probability proportional to their score
        parents = []

        # caclulate probability of selection normalized to 1
        scores = [x[1] for x in sorted_scores]
        scores = np.array(scores)
        scores = scores / np.sum(scores)
        print('Scores: \n', scores)
        # select compounds only
        compounds_only = [x[0] for x in sorted_scores]
        # in range of half of population size generate random number and select compound
        # Starting from the top of the population, keep adding the finesses to the partial sum P, till P<S
        # The individual for which P exceeds S is the chosen individual.

        # generate offsprings
        parents = np.random.choice(self.current_population, size=int(self.population_size / 2), p=scores, replace=False)
        parents_1 = parents[:len(parents) // 2]
        parents_2 = parents[len(parents) // 2:]
        offsprings_1 = [parent_1.cross(parent_2) for parent_1, parent_2 in zip(parents_1, parents_2)]
        offsprings_2 = [parent_2.cross(parent_1) for parent_1, parent_2 in zip(parents_1, parents_2)]
        offsprings = offsprings_1 + offsprings_2

        print(f'{len(offsprings)} offsprings generated from {len(parents)} parents.')
        print(f'Number unique parents: {len(set(parents))}')

        # in same fashion select subset of population to mutate
        mutate_rolls = np.random.rand(int(self.population_size / 2))
        # select best 20% of original population
        elite = sorted_scores[:int(self.population_size / 5)]
        elite = [x[0] for x in elite]
        # replace worst half of population with offsprings
        joined = list(parents) + offsprings + elite
        self.next_population = joined
        print('Next population size is: ', len(self.next_population))
        self._generate_population(self.next_population)
        self.history[self.current_iteration] = mean_scores
        return sorted_scores


GA = GeneticDocker(population_size=args.population_size)
GA.set_parts(pss, links, ligs)
for gen in range(args.generations):
    scores = GA.run_iteration()
    print(scores)
    print('-----------------------')
print(GA.history)
mean_scores_list = []
for iter in GA.history:
    print(iter)
    # calucate generation mean
    mean_scores = {compound: np.mean(scores) for compound, scores in GA.history[iter].items()}
    # sort by generation
    sorted_scores = sorted(mean_scores.items(), key=lambda x: x[1])
    mean_scores_list.append(np.mean([x[1] for x in sorted_scores]))
plt.plot(mean_scores_list)
plt.title('Mean score per generation')
plt.xlabel('Generation')
plt.ylabel('Mean score [kcal/mol]')
# plt.show()
# save fig with timestamp
plt.savefig(f'./{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}_evolution.png')

# plot best score per generation
best_scores_list = []
for iter in GA.history:
    print(iter)
    # calucate generation mean
    mean_scores = {compound: np.mean(scores) for compound, scores in GA.history[iter].items()}
    # sort by generation
    sorted_scores = sorted(mean_scores.items(), key=lambda x: x[1])
    best_scores_list.append(sorted_scores[0][1])
plt.plot(best_scores_list)
plt.title('Best score per generation')
plt.xlabel('Generation')
plt.ylabel('Best score [kcal/mol]')
plt.savefig(f'./{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}_best_score.png')
# plt.show()
# print best compound linage
# print_linage(sorted_scores[0][0])
# print('Best compound: ', sorted_scores[0][0].name)
# print('Best compound score: ', sorted_scores[0][1])
# if sorted_scores[0][0].parent_1 is not None and sorted_scores[0][0].parent_2 is not None:
#     print('Best compound parents: ', sorted_scores[0][0].parent_1.name, sorted_scores[0][0].parent_2.name)
#     print('Best compound parents scores: ', sorted_scores[0][0].parent_1.score, sorted_scores[0][0].parent_2.score)
#     if sorted_scores[0][0].parent_1.parent_1 != 'None':
#         print('Best compound grandparents: ', sorted_scores[0][0].parent_1.parent_1.name,
#               sorted_scores[0][0].parent_1.parent_2.name, sorted_scores[0][0].parent_2.parent_1.name,
#               sorted_scores[0][0].parent_2.parent_2.name)
#         print('Best compound grandparents scores: ', sorted_scores[0][0].parent_1.parent_1.score,
#               sorted_scores[0][0].parent_1.parent_2.score, sorted_scores[0][0].parent_2.parent_1.score,
#               sorted_scores[0][0].parent_2.parent_2.score)
# print('Is best compound mutated: ', sorted_scores[0][0].is_mutated)
# print('Is best compound offspring: ', sorted_scores[0][0].conjugate)
init_genes = GA.init_genes
init_genes['gen'] = 0
last_genes = GA.genes_overview(None)
last_genes['gen'] = args.generations - 1
merged = pd.concat([init_genes, last_genes])
merged.to_csv(f'./{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}_genes_overview.csv', index=False)
SCORES_DF.to_csv('chk.csv', index=True)
