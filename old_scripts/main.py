import argparse
import pandas as pd
import os, sys
import subprocess
from ..molecular_docking.smina_docking import docking_smina_single
from ..molecular_docking.gnina_docking import docking_gnina_single
# add molecular_docking to the path
sys.path.append('/home/anton/in_dev/Docking_tools')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("protein_file", help="protein file in pdbqt format")
    parser.add_argument("ligand_csv", help="ligand csv file")
    parser.add_argument("--num_modes", help="number of docking modes", default=5)
    parser.add_argument("--exhaustiveness", help="exhaustiveness of the docking", default=8)
    parser.add_argument("--autobox_ligand", help="ligand file in sdf format for autobox", default=None)
    parser.add_argument("--autobox_add", help="Distance to add to box in A", default=4)
    parser.add_argument("--seed", help="seed for random number generator", default=42)
    parser.add_argument("--box_def", help="box definition", default=None)
    args = parser.parse_args()
    print(args)

    # parse csv file
    df = pd.read_csv(args.ligand_csv)
    smiles = df['Conj_smiles'].values
    names = df['Molecule ChEMBL ID'].values
    for i in range(len(smiles)):
        print(names[i])
        print(i + 1)
        name = str(i + 1) + "-" + names[i]
        docking_smina_single(protein_file=args.protein_file,
                             ligand=smiles[i],
                             autobox_ligand=args.autobox_ligand,
                             exhaustiveness=args.exhaustiveness,
                             num_modes=args.num_modes,
                             docking_name=name,
                             out_dir="smina_out",
                             autobox_add=args.autobox_add,
                             box_def=args.box_def)


if __name__ == "__main__":
    main()
