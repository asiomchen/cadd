import os
import matplotlib.pyplot as plt
from utils.fixes import global_seed
import numpy as np
from create_conjugates.reaction import Reactor
import random
import json
import pickle

global_seed()

ANSWER = '''Vibrational spectroscopy is one of the most useful experimental
tools for study of hydrogen bonded complexes. So the information on calculated 
harmonic vibrational frequencies can be useful.'''
CHROMOSOME_LENGTH = len(ANSWER)
# genetic algorithm for simple string prediction
# available genes are lower and upper case letters, space and punctuation and special characters and polish letters
# fitness function is the number of correct characters in the string
GENES = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ 1234567890!@#$%^&*()_+-=[]{};:,./<>?|`~ąćęłńóśźżĄĆĘŁŃÓŚŹŻ'
MUTATION_RATE = 1 / CHROMOSOME_LENGTH


def random_chromosome():
    # create a random chromosome of length CHROMOSOME_LENGTH from GENES
    return "".join([random.choice(GENES) for _ in range(CHROMOSOME_LENGTH)])


def mutate(gene, mutation_rate=MUTATION_RATE):
    # mutate a single gene with a given mutation rate
    if random.random() < mutation_rate:
        return random.choice(GENES)
    return gene


class Individual:
    def __init__(self, string=None):
        self.string = string
        self.fitness = self._calc_fitness()
        self.age = 0

    def __repr__(self):
        return f'Individual(string={self.string}, fitness={self.fitness})'

    def _calc_fitness(self):
        fitness = len(ANSWER)
        for s, a in zip(self.string, ANSWER):
            if s == a:
                fitness -= 1
        return fitness ** 5 + fitness ** 3

    def crossover(self, other):
        '''
        Perform mating and produce new offspring
        '''

        # chromosome for offspring
        child_chromosome = []
        for gp1, gp2 in zip(self.string, other.string):

            # random probability
            prob = random.random()

            # if prob is less than 0.45, insert gene
            # from parent 1
            if prob < 0.45:
                child_chromosome.append(gp1)
            # if prob is between 0.45 and 0.90, insert
            # gene from parent 2
            elif prob < 0.90:
                child_chromosome.append(gp2)

            # otherwise insert random gene(mutate),
            # for maintaining diversity
            else:
                child_chromosome.append(mutate(gp1))

        # create new Individual(offspring) using
        # generated chromosome for offspring
        return Individual(string="".join(child_chromosome))


class GeneticAlgorythm:
    def __init__(self, population_size=100, mutation_rate=MUTATION_RATE):
        self.population_size = population_size
        self.mutation_rate = mutation_rate
        self.population = [Individual(random_chromosome()) for _ in range(population_size)]
        self.next_population = None
        self.generation = 0

    def run_iteration(self, elite_size=0.20, parents_size=0.5, strategy='default'):
        # run a single iteration of the genetic algorithm
        if self.next_population is not None:
            self.population = self.next_population
            self.next_population = None
        print(f'Generation {self.generation}')
        self.generation += 1
        # sort the population by fitness
        self.population.sort(key=lambda x: x.fitness)
        if self.population[0].fitness == 0:
            return self.population[0].fitness
        # add 1 to the age of all individuals
        for individual in self.population:
            individual.age += 1
        if len(self.population) > self.population_size:
            self.population = self.population[:self.population_size]

        fitness_sum = sum([ind.fitness for ind in self.population])
        # calculate the probability of each individual to be selected
        probabilities = [ind.fitness / fitness_sum for ind in self.population]
        # select the elite
        elite_size = int(self.population_size * elite_size)
        elite = self.population[:elite_size]
        score_rank = [1 / ind.fitness for ind in self.population]
        score_rank_sum = sum(score_rank)
        probabilities = [score_rank[i] / score_rank_sum for i in range(len(score_rank))]

        if strategy == 'default':
            parents = np.random.choice(self.population, size=int(self.population_size * parents_size), p=probabilities)
            parents_1 = parents[:len(parents) // 2]
            parents_2 = parents[len(parents) // 2:]
            offsprings_1 = [parent_1.crossover(parent_2) for parent_1, parent_2 in zip(parents_1, parents_2)]
            offsprings_2 = [parent_2.crossover(parent_1) for parent_1, parent_2 in zip(parents_1, parents_2)]
            self.next_population = list(elite) + list(offsprings_2) + list(parents) + list(offsprings_1)
        elif strategy == 'offspring':
            offspings = []
            while len(offspings) < self.population_size:
                parent_1 = np.random.choice(self.population, p=probabilities)
                parent_2 = np.random.choice(self.population, p=probabilities)
                offspings.append(parent_1.crossover(parent_2))
            self.next_population = list(elite) + list(offspings)

        return self.population[0].fitness


if __name__ == '__main__':
    print('Genetic Algorithm for Simple String Prediction, minimum fitness is {}'.format(CHROMOSOME_LENGTH))
    ind_1 = Individual(string=random_chromosome())
    ind_2 = Individual(string=random_chromosome())
    GA = GeneticAlgorythm(population_size=500, mutation_rate=MUTATION_RATE * 2)
    for _ in range(30000):
        fitness = GA.run_iteration(strategy='offspring')
        print(f'Best individual: {GA.population[0]}')
        print(f'Fitness: {fitness}')
        if fitness == 0:
            print('Found solution! Generation: {}'.format(GA.generation))
            break
