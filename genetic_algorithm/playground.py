import os
import matplotlib.pyplot as plt
from utils.fixes import global_seed
import numpy as np
import random

global_seed()


def calculate_rank_probabilities(population_size, a_max=1.2):
    a_min = 2 - a_max
    return [(a_max - (a_max - a_min) * (rank - 1) / (population_size - 1)) * 1 / population_size for rank in
            range(1, population_size + 1)]


def calculate_nfold(probabilities):
    return max(probabilities) / min(probabilities)


if __name__ == '__main__':
    print('Running simple test')
    test = calculate_nfold(calculate_rank_probabilities(50))
    print(test)
    a_mxs = np.linspace(1.1, 1.99, 30)
    folds = []
    for a_mx in a_mxs:
        print(a_mx)
        folds.append(calculate_nfold(calculate_rank_probabilities(50, a_mx)))
    plt.plot(a_mxs, folds)
    plt.ylim(0, 30)
    plt.show()
