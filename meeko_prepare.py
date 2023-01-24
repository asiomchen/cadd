from meeko import MoleculePreparation
from meeko import obutils
import argparse

parser = argparse.ArgumentParser(description='Prepare a molecule for docking')
parser.add_argument('--input', type=str, help='Input file')
parser.add_argument('--output', type=str, help='Output file', default='output.pdbqt')
# end

ars = parser.parse_args()

mol = obutils.load_molecule_from_file(ars.input)

preparator = MoleculePreparation(keep_nonpolar_hydrogens=False, macrocycle=True, hydrate=True)
preparator.prepare(mol)
preparator.show_setup()

output_pdbqt_file = ars.output
preparator.write_pdbqt_file(output_pdbqt_file)
