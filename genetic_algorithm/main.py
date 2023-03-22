import json
import os, sys

sys.path.insert(0, '..')
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
from utils.convert_alt import to_sdf
from utils.fixes import global_seed
from utils.extract_scores import extract_scores
from itertools import product
import numpy as np
from create_conjugates.reaction import Reactor
import random
import warnings
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
parser.add_argument('--parent_size', type=int, default=1, help='Parent size', required=False)
parser.add_argument('--mutation_rate', type=float, default=0.333, help='Mutation rate', required=False)
parser.add_argument('--selection', type=str, default='roulette', help='Selection method', required=False,
                    choices=['roulette', 'rank'])
parser.add_argument('--amax', type=float, default=1.2, help='a_max', required=False)
parser.add_argument('--duplicates', action='store_true', help='Allow duplicates', required=False)
parser.add_argument('--no-duplicates', action='store_false', help='Allow duplicates', required=False)
parser.set_defaults(duplicates=False)
parser.add_argument('--output', type=str, help='Output file', default='output.log')
parser.add_argument('--chkpt', type=str, help='Checkpoint file', default='./docked_mols', required=False)
parser.add_argument('--n_jobs', type=int, help='Number of jobs', default=2, required=False)
parser.add_argument('--seed', type=int, help='Random seed', default=42, required=False)
parser.add_argument('--plot_suffix', type=str, help='Suffix for plot', default='', required=False)
# end, parse arguments
args = parser.parse_args()
# load all compounds
data_path = '../data/parts/'
ligs = open(data_path + 'ligands_v4.smi').readlines()
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


def drop_duplicates(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


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
        return f'Compound({self.name}, score={self.score:.3f})'

    def __eq__(self, other):
        return (self.ps, self.link, self.lig) == (other.ps, other.link, other.lig)

    def __hash__(self):
        return hash((self.ps, self.link, self.lig))
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
        if os.path.exists(out_dir) is False:
            os.mkdir(out_dir)
        # check if molecule is already docked
        if os.path.exists(out_dir + '/' + self.name + '.pdbqt'):
            warnings.warn(f'Molecule {self.name} is already docked, check current algorithm')
        # docking command template using qvina-w
        cmd_template = "./qvina-w --receptor {} --ligand {} --num_modes 3 " \
                       "--exhaustiveness {} --seed 42 --out {} " \
                       "--center_x 42.8405 --center_y 31.0155 --center_z 32.3135 " \
                       "--size_x 34.751 --size_y 43.7889 --size_z 34.821"
        out_name = out_dir + '/' + self.name + '.pdbqt'
        # filling command template
        cmd = cmd_template.format(protein, ligand, ex, out_name)
        print(f'Entering docking for {self.name}...')
        # running docking without output
        return_code = subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if return_code != 0:
            print('Error in docking, smiles: ', self.conjugate)
            return self.score
        else:
            # if docking was successful, we extract scores from output pdbqt and mean is used as score
            scores = extract_scores(out_name)
            self.score = np.mean(scores)
            print(f'Exiting docking for {self.name}, score: {self.score}')
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
        sdf_name = to_sdf(self.conjugate, self.name, silent=True)
        if sdf_name is None:
            warnings.warn(f'Error in conversion to sdf for {self.name}')
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
    def __init__(self, population_size, selection='roulette', a_max=1.2):
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
        self.compounds_docked = set()
        self.selection = selection
        self.a_max = a_max
        if self.selection == 'rank':
            if self.a_max > 2 or self.a_max < 1:
                raise ValueError('a_max must be in range [1, 2]')
            self.rank_probabilities = self._calculate_rank_probabilities(self.population_size, a_max=self.a_max)
            print(f'Setting rank probabilities to:\n {self.rank_probabilities}')
            print(f'Sum of probabilities: {sum(self.rank_probabilities)}')
            print('Max probability: ', max(self.rank_probabilities))
            print('Min probability: ', min(self.rank_probabilities))
            if min(self.rank_probabilities) < 0.0001:
                print('WARNING: Min probability is less than 0.0001')
            else:
                diff = max(self.rank_probabilities) / min(self.rank_probabilities)
                print(f'Max / Min: {diff:.2f}')

    def _calculate_rank_probabilities(self, population_size, a_max=1.2):
        a_min = 2 - a_max
        return [(a_max - (a_max - a_min) * (rank - 1) / (population_size - 1)) * 1 / population_size for rank in
                range(1, population_size + 1)]

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

    def _init_population(self, size, unique=(not args.duplicates)):
        init_population = []
        names = []
        while len(init_population) != size:
            ps = random.sample(self.pss, 1)[0]
            link = random.sample(self.links, 1)[0]
            lig = random.sample(self.ligs, 1)[0]
            compound_name = f'{self.pss.index(ps)}_{self.links.index(link)}_{self.ligs.index(lig)}'
            if unique:
                if compound_name not in names:
                    init_population.append(Compound(ps, link, lig, name=compound_name))
                    names.append(compound_name)
            else:
                init_population.append(Compound(ps, link, lig, name=compound_name))
        print(f'Generating {len(init_population)} initial compounds...')
        self._generate_population(init_population)
        self.init_genes = self.genes_overview(init_population)
        return init_population

    def _generate_population(self, population):
        comp_count = 0
        to_generate = [comp for comp in population if comp.name + '.pdbqt' not in os.listdir('./generated_mols')]
        print(f'To generate before cleaning: {len(to_generate)}')
        # select only duplicates
        non_duplicates = set(to_generate)
        print(f'Number of duplicates: {len(to_generate) - len(non_duplicates)}')
        # out of duplicates, select only one instance
        to_generate = list(non_duplicates)
        comp_count = len(to_generate)
        with Pool(processes=10) as pool:
            pool.map(Compound.to_pdbqt, to_generate)
        print(f'Generated {comp_count} compounds  out of {len(population)}.')
        print('Generated compounds: ', [c.name for c in population])

    def set_parts(self, pss, links, ligs):
        self.pss = pss
        self.links = links
        self.ligs = ligs
        self.current_population = self._init_population(size=self.population_size)

    def run_iteration(self):
        global SCORES_DF
        mean_scores = {}
        mean_list = []
        print('Running iteration: ', self.current_iteration)

        if self.next_population is not None:
            self.current_population = self.next_population
            self.next_population = None
        self.compounds_docked = self.compounds_docked.union(set([c.name for c in self.current_population]))
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
                compound.score = SCORES_DF.loc[compound.name, 'score']
                mean_list.append(compound)
            if SCORES_DF is not None and compound.name not in SCORES_DF.index:
                docking_queue.append(compound)
        print(f'Queue length: {len(docking_queue)}')
        print(f'Queue: {docking_queue}')
        # run docking in parallel
        with Pool(args.n_jobs) as p:
            # use .dock() method to run docking
            scores = p.map(Compound.dock, docking_queue)
        print(f'Queue after: {docking_queue}')
        for compound, score in zip(docking_queue, scores):
            compound.score = score
            print(f'Compound: {compound.name}, score: {score}')
        print(f'Scores: {scores}')

        new_scores = {compound.name: compound.score for compound in docking_queue}
        list_only = mean_list + docking_queue

        sorted_scores = sorted(list_only, key=lambda x: x.score)
        for compound, score in zip(docking_queue, scores):
            mean_scores[compound] = score
        new_scores_df = pd.DataFrame.from_dict(new_scores, orient='index', columns=['score'])
        # add to SCORES_DF
        SCORES_DF = pd.concat([SCORES_DF, new_scores_df])
        # drop index duplicates
        SCORES_DF = SCORES_DF[~SCORES_DF.index.duplicated(keep='first')]

        if SCORES_DF is None:
            to_df = {c.name: score for c, score in mean_scores.items()}
            SCORES_DF = pd.DataFrame.from_dict(to_df, orient='index')
            SCORES_DF.columns = ['score']

        # apply elitism select only population_size best compounds
        print('Applying elitism...')
        self.current_population = sorted_scores[:self.population_size]
        print('Current population size: ', len(self.current_population))
        # select 3 mols with probability proportional to their score
        parents = []
        if self.selection == 'roulette':
            # caclulate probability of selection normalized to 1
            scores_ranking = [x.score for x in self.current_population]
            scores_ranking = np.array(scores_ranking)
            scores_ranking = scores_ranking / np.sum(scores_ranking)
            print('Scores: \n', scores_ranking)
        if self.selection == 'rank':
            scores_ranking = self.rank_probabilities
        # in range of half of population size generate random number and select compound
        # Starting from the top of the population, keep adding the finesses to the partial sum P, till P<S
        # The individual for which P exceeds S is the chosen individual.
        print('Selecting parents...')
        print(f'Population size: {len(self.current_population)}')
        print(f'Population: {self.current_population}')
        print(f'args.duplicates: {args.duplicates}')
        n_parents = int(self.population_size * args.parent_size)

        # generate offsprings
        if args.duplicates:
            print('Entering duplicates mode...')
            if n_parents % 2 != 0:
                n_parents += 1
            parents = np.random.choice(self.current_population, size=n_parents, p=scores_ranking,
                                       replace=True)
            parents_1 = parents[:len(parents) // 2]
            parents_2 = parents[len(parents) // 2:]
            offsprings_1 = [parent_1.cross(parent_2) for parent_1, parent_2 in zip(parents_1, parents_2)]
            offsprings_2 = [parent_2.cross(parent_1) for parent_1, parent_2 in zip(parents_1, parents_2)]
            offsprings = offsprings_1 + offsprings_2
            if len(offsprings) > n_parents:
                print('Too many offsprings, cutting off the rest.')
                offsprings = offsprings[:n_parents]

            print(f'{len(offsprings)} offsprings generated from {len(parents)} parents.')
            print(f'Number unique parents: {len(set(parents))} of {len(parents)}')
            print(f'Number unique offsprings: {len(set(offsprings))} of {len(offsprings)}')
            n_elite = int(self.population_size * args.elite_size)
            elite = sorted_scores[:n_elite]


        elif not args.duplicates:
            print('Entering no duplicates mode...')
            parents = np.random.choice(self.current_population, size=n_parents,
                                       p=scores_ranking, replace=True)
            parents_1 = parents[:len(parents) // 2]
            parents_2 = parents[len(parents) // 2:]
            offsprings_1 = [parent_1.cross(parent_2) for parent_1, parent_2 in zip(parents_1, parents_2)]
            offsprings_2 = [parent_2.cross(parent_1) for parent_1, parent_2 in zip(parents_1, parents_2)]
            initial_offsprings = drop_duplicates(offsprings_1 + offsprings_2)
            print(f'Initial offsprings: {len(initial_offsprings)} of {len(offsprings_1 + offsprings_2)}')
            while len(initial_offsprings) < int(self.population_size * args.parent_size):
                difference = n_parents - len(initial_offsprings)
                print(f'Adding {difference} offsprings...')
                if difference % 2 != 0:
                    difference += 1
                additional_parents = parents[:difference]
                print(f'Additional parents: {additional_parents}')
                additional_parents_1 = additional_parents[:len(additional_parents) // 2]
                additional_parents_2 = additional_parents[len(additional_parents) // 2:]
                additional_offsprings_1 = [parent_1.cross(parent_2) for parent_1, parent_2 in
                                           zip(additional_parents_1, additional_parents_2)]
                additional_offsprings_2 = [parent_2.cross(parent_1) for parent_1, parent_2 in
                                           zip(additional_parents_1, additional_parents_2)]
                additional_offsprings = additional_offsprings_1 + additional_offsprings_2
                initial_offsprings = drop_duplicates(initial_offsprings + additional_offsprings)
            offsprings = list(initial_offsprings)

            n_elite = int(self.population_size * args.elite_size)
            unique_elite = list(drop_duplicates(sorted_scores))
            elite = unique_elite[:n_elite]

        if len(offsprings) > n_parents:
            print('Too many offsprings, cutting off the rest.')
            offsprings = offsprings[:n_parents]
        print(f'{len(offsprings)} offsprings generated from {len(parents)} parents.')
        print(elite)
        # replace worst half of population with offsprings
        joined = offsprings + elite
        self.next_population = joined
        print('Next population size is: ', len(self.next_population))
        self._generate_population(self.next_population)
        self.history[self.current_iteration] = mean_scores
        self.current_iteration += 1
        return sorted_scores


# actual run part

GA = GeneticDocker(population_size=args.population_size, selection=args.selection, a_max=args.amax)
GA.set_parts(pss, links, ligs)
for gen in range(args.generations):
    scores = GA.run_iteration()
    print(scores)
    print('-----------------------')
print(f'Finishing algorithm after {args.generations} generations. {len(GA.compounds_docked)} compounds docked.')
mean_scores_list = []
generation_names = []
best_scores_list = []
best_5_scores_list = []
best_10_scores_list = []
best_20_scores_list = []
for iter in GA.history:
    print(iter)
    # calucate generation mean
    mean_scores = {compound: np.mean(scores) for compound, scores in GA.history[iter].items()}
    generation_compounds = [compound.name for compound in mean_scores.keys()]
    generation_names.append(' '.join(generation_compounds))
    # sort by generation
    sorted_scores = sorted(mean_scores.items(), key=lambda x: x[1])
    mean_scores_list.append(np.mean([x[1] for x in sorted_scores]))
    best_scores_list.append(sorted_scores[0][1])
    best_5_scores_list.append(np.mean([x[1] for x in sorted_scores[:5]]))
    best_10_scores_list.append(np.mean([x[1] for x in sorted_scores[:10]]))
    best_20_scores_list.append(np.mean([x[1] for x in sorted_scores[:20]]))
plt.plot(mean_scores_list)
# plot best score per generation
plt.plot(best_scores_list)
if args.generations >= 5:
    plt.plot(best_5_scores_list)
if args.generations >= 10:
    plt.plot(best_10_scores_list)
if args.generations >= 20:
    plt.plot(best_20_scores_list)
plt.title('Best score per generation')
plt.xlabel('Generation')
plt.ylabel('Best score [kcal/mol]')
plot_name = f'./GA_results/{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}_best_score_' + args.plot_suffix + '.png'
plt.legend(['Mean score', 'Best score', 'Best 5 score', 'Best 10 score', 'Best 20 score'])
plt.savefig(plot_name)
df_dict = {'gen': list(range(args.generations)), 'mean_score': mean_scores_list, 'compound': generation_names}
if args.generations >= 5:
    df_dict['best_5_score'] = best_5_scores_list
if args.generations >= 10:
    df_dict['best_10_score'] = best_10_scores_list
if args.generations >= 20:
    df_dict['best_20_score'] = best_20_scores_list
best_score = pd.DataFrame(df_dict)
best_score.to_csv(
    f'./GA_results/{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}_best_score_' + args.plot_suffix + '.csv', index=False)
# --------------------
init_genes = GA.init_genes
init_genes['gen'] = 0
last_genes = GA.genes_overview(None)
last_genes['gen'] = args.generations - 1
merged = pd.concat([init_genes, last_genes])
merged.to_csv(
    f'./GA_results/{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}_genes_overview_' + args.plot_suffix + '.csv',
    index=True)
SCORES_DF.to_csv('chk.csv', index=True)
with open(f'./GA_results/{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}_args_' + args.plot_suffix + '.json', 'w') as f:
    params = vars(args)
    params['docked_compounds'] = list(GA.compounds_docked)
    json.dump(params, f)
