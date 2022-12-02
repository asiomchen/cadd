import pandas as pd
import os
import subprocess

import pandas as pd

from utils.convert_alt import to_pdbqt, to_sdf


def docking_gnina_single(protein_file: str,
                         ligand: str,
                         autobox_ligand: str = None,
                         exhaustiveness: int = 8,
                         ligand_type: str = 'smiles',
                         num_modes: int = 1,
                         autobox_add: int = 4,
                         seed: int = 42,
                         docking_name: str = 'docked',
                         box_def: tuple = (
                                 11.338500022888184, -19.122999668121338, 181.1580047607422, 29.503000259399414,
                                 18.498000144958496,
                                 17.98199462890625),
                         out_dir: str = "gnina_output"):
    # searching for gnina in the tmp/gnina directory
    if not os.path.isfile("/tmp/gnina"):
        print("Gnina not found in the path")
        print('Downloading gnina...')
        cmd = 'wget https://github.com/gnina/gnina/releases/download/v1.0.1/gnina -P /tmp'
        subprocess.run(cmd, shell=True)
        cmd_2 = 'chmod +x /tmp/gnina'
        subprocess.run(cmd_2, shell=True)
    # creation of the output directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # ligand preparation
    if ligand_type == 'smiles':
        # out_name = convert_smiles.to_pqbqt(ligand, 'ligand', silent=False)
        out_name = to_pdbqt(ligand, docking_name)
    # ligand preparation end
    docking_name = os.path.join(out_dir, docking_name)
    # vina docking
    if autobox_ligand is not None:
        command = f"/tmp/gnina -r {protein_file} -l {out_name} " \
                  f"--autobox_ligand {autobox_ligand} --autobox_add {autobox_add} " \
                  f"--exhaustiveness {exhaustiveness} --num_modes {num_modes} --seed {seed} -o {docking_name}.sdf"
        print(command)
    else:
        command = f"/tmp/gnina -r {protein_file} -l {out_name} " \
                  f"--center_x {box_def[0]} --center_y {box_def[1]} --center_z {box_def[2]} " \
                  f"--size_x {box_def[3]} --size_y {box_def[4]} --size_z {box_def[5]} " \
                  f"--exhaustiveness {exhaustiveness} --num_modes {num_modes} --seed {seed} -o {docking_name}.sdf"
        print(command)

    gnina = subprocess.run(command, shell=True)
    # gnina docking end

    # removing the temporary file
    os.remove(out_name)


if __name__ == "__main__":
    os.chdir('/home/anton/in_dev/Docking_tools/master')
    df = pd.read_csv('./ligands/cox_chembl_cleaned_v3_conj.csv')
    smiles = df['Conj_smiles'].values
    names = df['Molecule ChEMBL ID'].values
    checkpoint = 0

    for i in range(checkpoint, len(smiles)):
        print(names[i])
        print(i + 1)
        name = str(i + 1) + "-" + names[i]
        docking_gnina_single(protein_file="./ligands/cox.pdbqt",
                             ligand=smiles[i],
                             autobox_ligand='./ligands/cox.pdbqt',
                             exhaustiveness=128,
                             ligand_type='smiles',
                             num_modes=5,
                             docking_name=name,
                             out_dir="gnina-128_output_conjugates")
