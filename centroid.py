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


def xyz_to_numpy(path):
    """
    Convert xyz file to numpy array
    Args:
        path: path to xyz file

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
            yield xyz


def centroid(xyz_generator):
    """
    Calculate centroid of a molecule
    Args:
        xyz_generator: generator of xyz files

    Returns: None

    """
    for xyz in xyz_generator:
        yield np.mean(xyz, axis=0)


if __name__ == '__main__':
    from_sdf_to_xyz('1-CHEMBL1077.sdf')
    print(xyz_to_numpy('1-CHEMBL1077.xyz'))
    for xyz in centroid(xyz_to_numpy('1-CHEMBL1077.xyz')):
        print(xyz)

    data = {'name': [],
            'pose': [],
            'centroid': []
            }
