import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools


class Reactor:
    # reaction def for basic conjugation, if linker is diamine
    # first converts amine linker and acid porphyrin to linked amine
    rcx_linker_to_ps_NN = '[CX4:0][NH2:1].[OH:2][C:3]=[0:4]>>[CX4:0][NH1:1][C:3]=[0:4]'
    # second reaction utilizes the same def, but used to conjugate linked porphyrin to ligand
    rcx_linked_add_ligand_NN = '[CX4:0][NH2:1].[OH:2][C:3]=[0:4]>>[CX4:0][NH1:1][C:3]=[0:4]'

    # for the sake of simplicity, this class assumes that ligands and porphyrins are carboxylic acids
    def __init__(self, linker_type='NN'):
        self.linker_type = linker_type

    def set_linker_type(self, linker_type):
        self.linker_type = linker_type

    def custom_reaction(self, reactants, rxn_smarts):
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        product = rxn.RunReactants(reactants)[0][0]
        Chem.SanitizeMol(product)
        return product

    def linker_to_ps(self, linker, ps):
        if self.linker_type == 'NN':
            rxn_smarts = self.rcx_linker_to_ps_NN
        return self.custom_reaction((linker, ps), rxn_smarts)

    def ligand_to_linked_ps(self, linked_ps, ligand):
        if self.linker_type == 'NN':
            rxn_smarts = self.rcx_linked_add_ligand_NN
        return self.custom_reaction((linked_ps, ligand), rxn_smarts)

    def conjugate(self, ps, linker, ligand):
        linked_ps = self.linker_to_ps(linker, ps)
        conjugate = self.ligand_to_linked_ps(linked_ps, ligand)
        return conjugate


if __name__ == '__main__':
    # ps = r"OC(C1=CC=C(/C2=C3C=CC(/C(C4=C(C=CC=C4F)F)=C5N/C(C=C\5)=C(C6=C(C=CC=C6F)F)\C7=N/C(C=C7)=C(C8=C(C=CC=C8F)F)\C9=CC=C2N9)=N\3)C=C1)=O"
    linker = 'NCCCCCCN'
    ligand = 'Cc1c(c2cc(ccc2n1C(=O)c3ccc(cc3)Cl)OC)CC(=O)O'
    # reactants = [Chem.MolFromSmiles(x) for x in (ps, linker, ligand)]
    #
    # product = Reactor().conjugate(*reactants)
    # print(Chem.MolToSmiles(product))
    # Chem.MolToSmiles(product, kekuleSmiles=False)
    # #save the molecule as a sdf
    # Chem.MolToMolFile(product, 'product_test.sdf')
    with open('/home/anton/PycharmProjects/cadd/data/parts/ps/ps.smi') as f:
        ps = f.read().splitlines()
        for i in ps:
            original_str = repr(i)[1:-1]

            reactants = [Chem.MolFromSmiles(x) for x in (original_str, linker, ligand)]
            # print(F'Reading: {Chem.MolToSmiles(ps)}')
            product = Reactor().conjugate(*reactants)
            print(F'Writing:  {Chem.MolToSmiles(product)}')
