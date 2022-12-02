import os
import subprocess
import warnings

import pandas as pd

from utils.convert_alt import to_pdbqt, to_sdf

warnings.filterwarnings("ignore")


def docking_smina_single(protein_file: str,
                         ligand: str,
                         autobox_ligand: str,
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
                         out_dir: str = "vina_output"):
    # searching for smina in the path
    if not os.path.isfile("smina.static"):
        print("Smina not found in the path. Downloading...")
        subprocess.call(["wget", "https://sourceforge.net/projects/smina/files/smina.static"])
        # make smina executable
        subprocess.call(["chmod", "+x", "smina.static"])
        print("Smina downloaded.")
    # creation of the output directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # ligand preparation
    if ligand_type == 'smiles':
        out_name = to_sdf(ligand, 'ligand')
        fix_uppercase(out_name)
    if ligand_type == 'pdbqt':
        out_name = ligand
    if ligand_type == 'sdf':
        out_name = ligand

    # ligand preparation end
    docking_name = os.path.join(out_dir, docking_name)
    # smina docking
    if autobox_ligand is not None:
        command = f"./smina.static -r {protein_file} -l {out_name} " \
                  f"--autobox_ligand {autobox_ligand} --autobox_add {autobox_add} --num_modes {num_modes}  " \
                  f"--exhaustiveness {exhaustiveness} --seed {seed} -o {docking_name}.sdf"
        print(command)
    else:
        command = f"./smina.static -r {protein_file} -l {out_name} " \
                  f"--center_x {box_def[0]} --center_y {box_def[1]} --center_z {box_def[2]} " \
                  f"--size_x {box_def[3]} --size_y {box_def[4]} --size_z {box_def[5]} --num_modes {num_modes}  " \
                  f"--exhaustiveness {exhaustiveness} --seed {seed} -o {docking_name}.sdf"
        print(command)
        return 1

    smina = subprocess.run(command, shell=True)
    # smina docking end

    # removing the temporary file
    os.remove(out_name)


def fix_uppercase(file, pattern=[' Cl', ' Br']):
    """
    In case of the ligand pdbqt file has second uppercase letter in the atom names, which is unapropriate for vina
    Pattern ['CL', 'BR'] will be replaced by ['Cl', 'Br']
    Pattern always occurs at the end of the line.
    """
    with open(file, 'r') as f:
        lines = f.readlines()
    with open(file, 'w') as f:
        for line in lines:
            for p in pattern:
                line = line.replace(p.upper(), p)
            f.write(line)


if __name__ == "__main__":
    os.chdir("/home/anton/in_dev/Docking_tools/master")
    df = pd.read_csv('/master/ligands/cox_chembl_cleaned_v3_conj.csv')
    smiles = df['Conj_smiles'].values
    names = df['Molecule ChEMBL ID'].values
    checkpoint = 0

    for i in range(checkpoint, len(smiles)):
        print(names[i])
        print(i + 1)
        name = str(i + 1) + "-" + names[i]
        docking_smina_single(protein_file="cox.pdbqt",
                             ligand=smiles[i],
                             autobox_ligand='cox.pdbqt',
                             exhaustiveness=8,
                             ligand_type='smiles',
                             num_modes=5,
                             docking_name=name,
                             out_dir="smina_conjugates", )
