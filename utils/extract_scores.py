import pandas as pd
from utils.centroid import *
from typing import Union, List, Generator, Dict


class Extractor:
    '''
    class for extracting scores from pdbqt files
    '''

    def __init__(self, mode='folder', format='pdbqt', data_path=None):
        self.extracted_data = None
        self.mode = mode
        self.data_path = data_path

    def set_format(self, format):
        self.format = format

    def set_data_path(self, data_path):
        self.data_path = data_path

    def extract(self):
        if self.mode == 'folder':
            self.extracted_data = self._extract_folder()
            return self.extracted_data
        elif self.mode == 'file':
            self.extracted_data = self._extract_file()
            return self.extracted_data
        else:
            raise ValueError('mode should be either "folder" or "file"')

    def _extract_folder(self):
        if self.data_path is None:
            raise ValueError('data_path is not set')
        else:
            data = extract_scores_from_dir(self.data_path)
            return data

    def _extract_file(self):
        if self.data_path is None:
            raise ValueError('data_path is not set')
        else:
            data = extract_scores_from_file(self.data_path)
            return data


class Point:
    '''
    base class for representing a point in 3D space
    '''

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return f'Point({self.x}, {self.y}, {self.z})'


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
        yield Point(xyz[0], xyz[1], xyz[2])


def extract_scores_from_file(path: str) -> pd.DataFrame:
    '''
    Extracts the name, scores and coordinates of the centroid from a file
    :param path: str
    :return: dict, which could be converted to a pandas DataFrame
    '''
    out = {'name': [], 'pose': [], 'score': [], 'centroid': []}
    scores = extract_scores(path)
    centroids = extract_coords(path)
    name = os.path.basename(path).split('.')[0]
    for i, score in enumerate(scores):
        out['name'].append(name)
        out['pose'].append(i)
        out['score'].append(score)
        out['centroid'].append(next(centroids))
    return pd.DataFrame(out)


def extract_scores_from_dir(path, type='pdbqt') -> Dict:
    '''
    Extracts the name, scores and coordinates of the centroid from a directory of docked files
    :param path: str
    :param type: str (default: pdbqt)
    :return: dict, which could be converted to a pandas DataFrame
    '''
    out = None
    for file in os.listdir(path):
        if out is None:
            out = extract_scores_from_file(os.path.join(path, file))
        else:
            out = pd.concat([out, extract_scores_from_file(os.path.join(path, file))])
    return out



if __name__ == '__main__':
    path = '../results/ares_qvina/docked_pocket_16/'
    test = '/home/anton/PycharmProjects/cadd/results/ares_qvina/docked_pocket_16/1_1_925.pdbqt'
    # scores = extract_scores_from_dir(path)
    # print(pd.DataFrame(scores))
    # xyz_path = from_pdbqt_to_xyz(test)
    # for xyz in centroid(xyz_to_numpy(xyz_path)):
    #     print(xyz)
    e = Extractor(mode='file', data_path=test)
    e2 = Extractor(mode='folder', data_path=path)
    print(e.extract().head())
    print(e2.extract().head())
