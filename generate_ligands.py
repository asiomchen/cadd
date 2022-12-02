import argparse

from utils.convert_alt import to_sdf
import pandas as pd
import os


def main():
    '''
    Converts a csv file with smiles to sdf files folder for docking
    :return:
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-ligand_csv", help="ligand csv file")
    parser.add_argument("-smiles_c", help="smiles column name", default="Conj_smiles")
    parser.add_argument("-out_dir", help="output directory", default="ligands")
    args = parser.parse_args()
    print(args)

    # parse csv file
    df = pd.read_csv(args.ligand_csv)
    smiles = df[args.smiles_c].values
    names = df['Molecule ChEMBL ID'].values

    # create output directory
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    for i in range(len(smiles)):
        print(names[i])
        print(i + 1)
        name = str(i + 1) + "-" + names[i]
        output_path = os.path.join(args.out_dir, name)
        to_sdf(smiles[i], output_path)


if __name__ == "__main__":
    main()
