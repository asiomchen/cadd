import os
import subprocess

import numpy as np

def from_sdf_to_xyz(path, silent: bool = True):
    """
    Convert sdf file to xyz file
    Args:
        path: path to sdf file

    Returns: None

    """
    path_to_xyz = path.replace('.sdf', '.xyz')
    cmd = "obabel " + path + " -O " + path_to_xyz
    if not silent:
        print(cmd)
    obabel_return_code = subprocess.run(cmd, shell=True).returncode
    return path_to_xyz


def from_pdbqt_to_xyz(path, silent: bool = True):
    """
    Convert pdbqt file to xyz file
    Args:
        path: path to sdf file

    Returns: None

    """
    path_to_xyz = os.path.join(os.path.dirname(path), os.path.basename(path).replace('.pdbqt', '.xyz'))
    cmd = "obabel " + path + " -O " + path_to_xyz
    if not silent:
        print(cmd)
    obabel_return_code = subprocess.run(cmd, shell=True).returncode
    return path_to_xyz


def xyz_to_numpy(path, remove_file: bool = True):
    """
    Convert xyz file to numpy array
    Args:
        path: path to xyz file
        remove_file: bool (default: True)

    Returns: None

    """
    with open(path, 'r') as f:
        lines = f.readlines()
        n_atoms = int(lines[0])
        xyz = np.zeros((n_atoms, 3))
        n_mols = len(lines) // (n_atoms + 2)
        print(f'Number of molecules in xyz file: {n_mols}')
        for i in range(n_mols):
            for j in range(n_atoms):
                xyz[j] = lines[i * (n_atoms + 2) + 2 + j].split()[1:4]
            yield xyz.copy()


def centroid(xyz_generator):
    """
    Calculate centroid and coordinates of the molecule
    Args:
        xyz_generator: generator of xyz files

    Returns: None

    """
    centroids = []
    coords = []
    for xyz in xyz_generator:
        centroid = np.mean(xyz, axis=0)
        centroids.append(centroid)
        coords.append(xyz)
    return centroids, coords


def get_box(xyz, extend=2):
    """
    Get box around coordinates
    Args:
        xyz: numpy array of xyz coordinates
        extend: int (default: 2)

    Returns: None

    """
    min_x = np.min(xyz[:, 0]) - extend
    max_x = np.max(xyz[:, 0]) + extend
    min_y = np.min(xyz[:, 1]) - extend
    max_y = np.max(xyz[:, 1]) + extend
    min_z = np.min(xyz[:, 2]) - extend
    max_z = np.max(xyz[:, 2]) + extend

    size_x = max_x - min_x
    size_y = max_y - min_y
    size_z = max_z - min_z
    center_x = (max_x + min_x) / 2
    center_y = (max_y + min_y) / 2
    center_z = (max_z + min_z) / 2

    return [center_x, center_y, center_z], [size_x, size_y, size_z]


if __name__ == '__main__':
    os.chdir('../data/test/')
    file = '1_1_359.pdbqt'
    conv = from_pdbqt_to_xyz(file)
    xyz = xyz_to_numpy(conv)
    res = centroid(xyz)
    print(res[0][0])

    # center, size = get_box(xyz)
    # print('Box:', [*center, *size])
    # all_boxes.append(['gridbox', *center, *size])

    # for xyz in centroid(xyz_to_numpy(file.replace('.pdbqt', '.xyz'))):
    #     print(xyz)
    #
    # data = {'name': [],
    #         'pose': [],
    #         'centroid': []
    #         }
