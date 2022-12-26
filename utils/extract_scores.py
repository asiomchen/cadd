import pandas as pd
from utils.centroid import *
from typing import Union, List, Generator, Dict


def extract_scores(path: str) -> List[float]:
    '''
    Extracts the scores from a pdbqt file
    :param path: str
    :return: list
    '''
    with open(path, 'r') as f:
        lines = f.readlines()
    scores = []
    for line in lines:
        if 'REMARK VINA RESULT:' in line:
            scores.append(float(line.split()[-3]))
    return scores


def extract_coords(path: str) -> Generator:
    '''
    Extracts the coordinates of the centroid from a pdbqt file
    :param path: str
    :return: Generator
    '''
    xyz_path = from_pdbqt_to_xyz(path)
    for xyz in centroid(xyz_to_numpy(xyz_path)):
        yield xyz


def extract_scores_from_dir(path, type='pdbqt') -> Dict:
    '''
    Extracts the name, scores and coordinates of the centroid from a directory of docked files
    :param path: str
    :param type: str (default: pdbqt)
    :return: dict, which could be converted to a pandas DataFrame
    '''

    out = {'name': [], 'pose': [], 'score': [], 'centroid': []}
    for file in os.listdir(path):
        if file.endswith('.pdbqt'):
            scores = extract_scores(os.path.join(path, file))
            centroids = extract_coords(os.path.join(path, file))
            for i, score in enumerate(scores):
                out['name'].append(file.split('.')[0])
                out['pose'].append(i)
                out['score'].append(score)
                out['centroid'].append(next(centroids))

    return out


if __name__ == '__main__':
    path = '/results/ares_qvina/docked_pocket_16/'
    test = '/home/anton/PycharmProjects/cadd/results/ares_qvina/docked_pocket_16/1_1_925.pdbqt'
    scores = extract_scores_from_dir(path)
    print(pd.DataFrame(scores))
    xyz_path = from_pdbqt_to_xyz(test)
    for xyz in centroid(xyz_to_numpy(xyz_path)):
        print(xyz)
