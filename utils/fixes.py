import random
import os
import numpy as np


def global_seed(seed=42):
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    np.random.seed(seed)
    print("Global seed set to {}".format(seed))
