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

parser = argparse.ArgumentParser()
# add arguments: population size, number of generations, elite size, parent size
parser.add_argument('--population_size', type=int, default=50, help='Population size', required=False)
parser.add_argument('--generations', type=int, default=20, help='Number of generations', required=False)
parser.add_argument('--elite_size', type=int, default=0.2, help='Elite size', required=False)
parser.add_argument('--parent_size', type=int, default=0.5, help='Parent size', required=False)
parser.add_argument('--mutation_rate', type=float, default=0.333, help='Mutation rate', required=False)
parser.add_argument('--output', type=str, help='Output file', default='output.log')
parser.add_argument('--chkpt', type=str, help='Checkpoint file', default='./docked_mols', required=False)
parser.add_argument('--n_jobs', type=int, help='Number of jobs', default=2, required=False)
parser.add_argument('--seed', type=int, help='Random seed', default=42, required=False)
parser.add_argument('--plot_suffix', type=str, help='Suffix for plot', default='', required=False)
# end, parse arguments
args = parser.parse_args()
# load all compounds
data_path = '../data/parts/'
ligs = open(data_path + 'ligands_v2.smi').readlines()
links = open(data_path + 'linker.smi').readlines()
pss = open(data_path + 'ps.smi').readlines()
# number of compounds possible to create
all_compounds = list(product(ligs, links, pss))
print(len(all_compounds))
# checkpoint is either a csv file or a folder with docked molecules
SCORES_DF = None
if args.chkpt and 'csv' in args.chkpt:
    SCORES_DF = pd.read_csv(args.chkpt, index_col=0)
    print(f'Loaded {len(SCORES_DF)} molecules from {args.chkpt}')
elif args.chkpt and './' in args.chkpt:
    SCORES_DF = retrive_chkpt(args.chkpt)
    print(f'Loaded {len(SCORES_DF)} molecules from folder {args.chkpt}')

# set random seed
global_seed(args.seed)


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


# ii theory all compound are cleaned, but to be sure we sanitize them
ligs = list(map(lambda x: global_sanitize(x), ligs))
links = list(map(lambda x: global_sanitize(x), links))
pss = list(map(lambda x: global_sanitize(x), pss))


class Compound:
    '''
    Base class for fragment-base compound descriptor.
    '''

    def __init__(self, ps, link, lig, name='compound', parent_1=None, parent_2=None, is_mutated=False):
        self.ps = ps
        self.link = link
        self.lig = lig
        self._conjugate = None
        self.name = name
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

    def dock(self, protein='../data/prots/cox2.pdbqt', ex=32, centroid=(42.84, 31.02, 32.31, 34, 75, 43.79, 34.82),
             out_dir='./docked_mols'):
        # now docking assuming that we only COX-2 and best pocket is used (to be changed)
        # if any problem with conjugate generation or/and docking, we return 0 as neutral score
        self.score = 0
        if self.conjugate is None:
            return self.score
        ligand = './generated_mols/' + self.name + '.pdbqt'
        if ligand is None:
            print('Error in generation, smiles: ', self.conjugate)
            return self.score
        else:
            ligand = './generated_mols/' + self.name + '.pdbqt'
        # docking command template using qvina-w
        cmd_template = "./qvina-w --receptor {} --ligand {} --num_modes 3 " \
                       "--exhaustiveness {} --seed 42 --out {} " \
                       "--center_x 42.8405 --center_y 31.0155 --center_z 32.3135 " \
                       "--size_x 34.751 --size_y 43.7889 --size_z 34.821"
        out_name = out_dir + '/' + self.name + '.pdbqt'
        # filling command template
        cmd = cmd_template.format(protein, ligand, ex, out_name)
        print(f'Entering docking for {self.name}')
        # running docking without output
        return_code = subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if return_code != 0:
            print('Error in docking, smiles: ', self.conjugate)
            return self.score
        else:
            # if docking was successful, we extract scores from output pdbqt and mean is used as score
            scores = extract_scores(out_name)
            print(f'Exiting docking for {self.name}, scores: {scores}')
            self.score = np.mean(scores)
            return self.score

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
        return offspring

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
            try:
                os.remove(sdf_name)
            except UserWarning:
                print(f'Something went wrong with {sdf_name} removal')
                if os.path.exists(sdf_name):
                    print(f'{sdf_name} still exists')
                if os.path.exists(self.name + '.pdbqt'):
                    print(f'{self.name}.pdbqt was created successfully')
                else:
                    print(f'{self.name}.pdbqt was not created')
        os.chdir('../')
        # return None if conversion was not successful
        if return_code != 0:
            print('Error in meeko, smiles: ', self.conjugate)
            return None
        else:
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
            if SCORES_DF is None:
                print('SCORES_DF is None')
            if SCORES_DF is not None and compound.name in SCORES_DF.index:
                mean_scores[compound] = SCORES_DF.loc[compound.name, 'score']
            if SCORES_DF is not None and compound.name not in SCORES_DF.index:
                docking_queue.append(compound)
        print(f'Queue length: {len(docking_queue)}')
        print(f'Queue: {docking_queue}')
        # run docking in parallel
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
        print(f'Number unique parents: {len(set([parent.name for parent in parents]))} of {len(parents)}')
        print(
            f'Number unique offsprings: {len(set([offspring.name for offspring in offsprings]))} of {len(offsprings)}')

        # in same fashion select subset of population to mutate
        # select best 20% of original population
        elite = sorted_scores[:int(self.population_size * args.elite_size)]
        elite = [x[0] for x in elite]
        # replace worst half of population with offsprings
        joined = list(parents) + offsprings + elite
        self.next_population = joined
        print('Next population size is: ', len(self.next_population))
        self._generate_population(self.next_population)
        self.history[self.current_iteration] = mean_scores
        self.current_iteration += 1
        return sorted_scores


# actual run part

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
plot_name = f'./{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}_evolution_' + args.plot_suffix + '.png'
plt.savefig(plot_name)
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
plot_name = f'./{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}_best_score_' + args.plot_suffix + '.png'
plt.savefig(plot_name)
# --------------------
init_genes = GA.init_genes
init_genes['gen'] = 0
last_genes = GA.genes_overview(None)
last_genes['gen'] = args.generations - 1
merged = pd.concat([init_genes, last_genes])
merged.to_csv(f'./{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}_genes_overview.csv', index=False)
SCORES_DF.to_csv('chk.csv', index=True)
