import os
import matplotlib.pyplot as plt
from utils.convert_alt import to_pdbqt, to_sdf
from utils.fixes import global_seed
from utils.extract_scores import extract_scores
from itertools import product
import numpy as np
from create_conjugates.reaction import Reactor
import random
from rdkit import Chem
import subprocess
import json
import pickle

global_seed()

data_path = '../data/parts/'
ligs = open(data_path + 'ligands_v2.smi').readlines()
links = open(data_path + 'linker.smi').readlines()
pss = open(data_path + 'ps.smi').readlines()
all_compounds = list(product(ligs, links, pss))
print(len(all_compounds))


def global_sanitize(smiles):
    smiles = smiles.split('.')[0]
    smiles = smiles.replace('[2H]', '[H]')
    smiles = smiles.replace('[3H]', '[H]')
    return smiles


ligs = list(map(lambda x: global_sanitize(x), ligs))
links = list(map(lambda x: global_sanitize(x), links))
pss = list(map(lambda x: global_sanitize(x), pss))


def MolWithoutIsotopesToSmiles(mol):
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
    Base class for fragment-base compound descriptor
    '''

    def __init__(self, ps, link, lig, name='compound', parent_1=None, parent_2=None, is_mutated=False):
        self.ps = self._sanitize(ps)
        self.link = self._sanitize(link)
        self.lig = self._sanitize(lig)
        self.conjugate = self._reaction()
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

    def __str__(self):
        return f'Compound(name={self.name}, parent_1={self.parent_1.name if self.parent_1 is not None else None}, parent_2={self.parent_2.name if self.parent_2 is not None else None}, is_mutated={self.is_mutated})'

    def __repr__(self):
        return f'Compound({self.name})'

    def __eq__(self, other):
        return self.ps == other.ps and self.link == other.link and self.lig == other.lig

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
        self.conjugate = self._reaction()

    def dock(self, protein='../data/prots/cox2.pdbqt', ex=32, centroid=(42.84, 31.02, 32.31, 34, 75, 43.79, 34.82),
             out_dir='./docked_mols'):
        scores = []
        if self.conjugate is None:
            self.score = 0
            return [0, 0, 0]
        if self.is_docked:
            scores = extract_scores(out_dir + '/' + self.name + '.pdbqt')
            self.score = np.mean(scores)
            return scores

        if self.is_generated:
            ligand = './generated_mols/' + self.name + '.pdbqt'
        else:
            ligand = self.to_pdbqt()
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
            return [0, 0, 0]
        else:
            scores = extract_scores(out_name)
            print(f'Exiting docking for {self.name}, scores: {scores}')
            return scores

    def _sanitize(self, smiles):
        smiles = smiles.split('.')[0]
        smiles = smiles.replace('[2H]', '[H]')
        smiles = smiles.replace('[3H]', '[H]')
        return smiles

    def _reaction(self):
        smiles_compounds = [self.ps, self.link, self.lig]
        mols = list(map(Chem.MolFromSmiles, smiles_compounds))
        reactor = Reactor()
        mol = reactor.conjugate(*mols)
        if mol is None:
            return None
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol)

    def cross(self, other, mutation_rate=0.333):
        '''
        Cross two compounds
        :param other:
        :return:
        '''
        p = [0.5, 0.5]
        # generate mutation mask
        mask = np.random.choice([0, 1], size=3, p=[1 - mutation_rate, mutation_rate])
        ps = np.random.choice([self.ps, other.ps], p=p)
        is_mutated = any(mask)
        if mask[0] == 1:
            ps = np.random.choice(pss)
        link = np.random.choice([self.link, other.link], p=p)
        if mask[1] == 1:
            link = np.random.choice(links)
        lig = np.random.choice([self.lig, other.lig], p=p)
        if mask[2] == 1:
            lig = np.random.choice(ligs)
        name = f'{pss.index(ps)}_{links.index(link)}_{ligs.index(lig)}'
        offspring = Compound(ps, link, lig, name=name, parent_1=self, parent_2=other, is_mutated=is_mutated)
        if offspring.conjugate is None:
            print('Error in reaction for compound in .cross(): ', name)
            print('Scoring will be skipped')
            offspring.name = 'ERROR_' + name
        offspring.reset()
        return offspring

    def mutate(self, mutation_rate=0.1):
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
        # convert smiles to pdbqt
        if self.conjugate is None:
            return None
        os.chdir(path)
        sdf_name = to_sdf(self.conjugate, self.name)
        meeko_cmd = f'mk_prepare_ligand.py -i {sdf_name} -o {self.name}.pdbqt'
        return_code = subprocess.call(meeko_cmd, shell=True)
        # remove sdf file
        os.remove(sdf_name)
        os.chdir('../')
        if return_code != 0:
            print('Error in meeko, smiles: ', self.conjugate)
            return None
        else:
            self.is_generated = True
            return self.name + '.pdbqt'


class GeneticDocker:
    def __init__(self, population_size):
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

    def _update_checkpoint(self, filename='checkpoint.csv'):
        if self.ligs is None or self.links is None or self.pss is None:
            raise ValueError('Parts are not set')
        # check if file exists if not create it if yes read it
        if not os.path.exists(filename):
            with open(filename, 'w') as f:
                f.write('0')
        else:
            pass

    def _init_population(self, size):
        init_population = []
        for i in range(size):
            ps = random.sample(self.pss, 1)[0]
            link = random.sample(self.links, 1)[0]
            lig = random.sample(self.ligs, 1)[0]
            compound_name = f'{self.pss.index(ps)}_{self.links.index(link)}_{self.ligs.index(lig)}'

            init_population.append(Compound(ps, link, lig, name=compound_name))
        self._generate_population(init_population)
        return init_population

    def _generate_population(self, population):
        comp_count = 0
        for compound in population:
            if compound.name + '.pdbqt' not in os.listdir('./generated_mols'):
                compound.to_pdbqt()
                comp_count += 1
        print(
            f'Generated {comp_count} compounds  out of {len(population)}. {len(population) - comp_count} already exist')

    def run(self, protein, ex=32, centroid=(42.84, 31.02, 32.31, 34, 75, 43.79, 34.82), out_dir='out'):
        scores = []

        cmd_template = "./qvina-w --receptor {} --ligand {} --num_modes $num_modes " \
                       "--exhaustiveness {} --seed 42 --out {} " \
                       "--center_x 42.8405 --center_y 31.0155 --center_z 32.3135 " \
                       "--size_x 34.751 --size_y 43.7889 --size_z 34.821"
        for compound in self.current_population:
            out_name = out_dir + '/' + compound.name
            cmd = cmd_template.format(protein, ex, out_name)

    def set_parts(self, pss, links, ligs):
        self.pss = pss
        self.links = links
        self.ligs = ligs
        self.current_population = self._init_population(size=self.population_size)

    def run_iteration(self):
        all_scores = {}
        print('Running iteration: ', self.current_iteration)

        if self.next_population is not None:
            self.current_population = self.next_population
            self.next_population = None
        print('Current compound names are: ', [c.name for c in self.current_population])
        print('{} compounds to dock'.format(sum([not c.is_docked for c in self.current_population])))
        for compound in self.current_population:
            if not compound.is_docked:
                compound.reset()
            # print('Docking compound: ', compound.name)
            # print('Compound is docked: ', compound.is_docked)
            all_scores[compound] = compound.dock()
        # apply elitism select only population_size best compounds

        self.current_iteration += 1
        mean_scores = {compound: np.mean(scores) for compound, scores in all_scores.items()}
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
        # select compounds only
        compounds_only = [x[0] for x in sorted_scores]
        # in range of half of population size generate random number and select compound
        # Starting from the top of the population, keep adding the finesses to the partial sum P, till P<S
        # The individual for which P exceeds S is the chosen individual.
        all_parents = np.random.rand(int(self.population_size))
        # split parents into 2 groups
        parents_1 = all_parents[:int(self.population_size / 2)]
        parents_2 = all_parents[int(self.population_size / 2):]

        # generate offsprings
        offsprings = []
        parent_compounds = []
        for i in range(int(self.population_size / 2)):
            parent_1 = compounds_only[np.where(np.cumsum(scores) > parents_1[i])[0][0]]
            parent_2 = compounds_only[np.where(np.cumsum(scores) > parents_2[i])[0][0]]
            offsprings.append(parent_1.cross(parent_2))
            parent_compounds.append(parent_1)
            parent_compounds.append(parent_2)
        print(f'{len(offsprings)} offsprings generated from {len(all_parents)} parents')

        # in same fashion select subset of population to mutate
        mutate_rolls = np.random.rand(int(self.population_size / 2))
        pre_mutated = parents.copy()

        mutated = [c.mutate(mutation_rate=0.2) for c in pre_mutated]
        # select best 20% of original population
        elite = sorted_scores[:int(self.population_size / 5)]
        elite = [x[0] for x in elite]
        # replace worst half of population with offsprings
        joined = parent_compounds + offsprings + elite
        # while len(joined) < self.population_size:
        #     joined.append(np.random.choice([x[0] for x in sorted_scores], p=scores))
        #     joined = [x for x in joined if x.name.contains('ERROR') is False]
        self.next_population = joined
        print('Next population size is: ', len(self.next_population))
        self._generate_population(self.next_population)
        self.history[self.current_iteration] = all_scores
        return sorted_scores


# comp = Compound(pss[0], links[0], ligs[0], name='0_0_0')
# comp_2 = Compound(pss[5], links[3], ligs[100], name='5_3_100')
# comp_3 = comp.cross(comp_2)
# print(comp_3)

def print_linage(compound):
    if compound.parent_1 is not None:
        print_linage(compound.parent_1)
    if compound.parent_2 is not None:
        print_linage(compound.parent_2)
    else:
        print(compound.name)
        print('-------')


GA = GeneticDocker(population_size=50)
GA.set_parts(pss, links, ligs)
for i in range(10):
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
plt.show()
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
plt.show()
# print best compound linage
# print_linage(sorted_scores[0][0])
print('Best compound: ', sorted_scores[0][0].name)
print('Best compound score: ', sorted_scores[0][1])
if sorted_scores[0][0].parent_1 is not None and sorted_scores[0][0].parent_2 is not None:
    print('Best compound parents: ', sorted_scores[0][0].parent_1.name, sorted_scores[0][0].parent_2.name)
    print('Best compound parents scores: ', sorted_scores[0][0].parent_1.score, sorted_scores[0][0].parent_2.score)
    if sorted_scores[0][0].parent_1.parent_1 != 'None':
        print('Best compound grandparents: ', sorted_scores[0][0].parent_1.parent_1.name,
              sorted_scores[0][0].parent_1.parent_2.name, sorted_scores[0][0].parent_2.parent_1.name,
              sorted_scores[0][0].parent_2.parent_2.name)
        print('Best compound grandparents scores: ', sorted_scores[0][0].parent_1.parent_1.score,
              sorted_scores[0][0].parent_1.parent_2.score, sorted_scores[0][0].parent_2.parent_1.score,
              sorted_scores[0][0].parent_2.parent_2.score)
print('Is best compound mutated: ', sorted_scores[0][0].is_mutated)
print('Is best compound offspring: ', sorted_scores[0][0].conjugate)
