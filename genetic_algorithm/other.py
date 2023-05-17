import numpy as np
import pandas as pd
import os, sys

sys.path.insert(0, '..')
from utils.extract_scores import extract_scores


def retrive_chkpt(path):
    chk = {}
    for file in os.listdir(path):
        if file.endswith('.pdbqt'):
            base_name = file.split('.')[0]
            chk[base_name] = np.mean(extract_scores(path + '/' + file))

    return pd.DataFrame.from_dict(chk, orient='index', columns=['score'])


if __name__ == '__main__':
    path = sys.argv[1]
    chk = retrive_chkpt(path)
    chk.to_csv('chk.csv')
    print(chk)
