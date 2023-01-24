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

    rcx_linker_to_ps_OO = '[CX4:0][OH:1].[OH:2][C:3]=[0:4]>>[CX4:0][0:1][C:3]=[0:4]'
    rcx_linked_add_ligand_OO = '[CX4:0][OH:1].[OH:2][C:3]=[0:4]>>[CX4:0][0:1][C:3]=[0:4]'

    # for the sake of simplicity, this class assumes that ligands and porphyrins are carboxylic acids
    def __init__(self, linker_type='NN', mode='Mol'):
        self.linker_type = None
        self.mode = mode

    def set_linker_type(self, linker_type):
        '''
        This function sets the linker type for the reaction.
        :param linker_type:
        :return:
        '''
        self.linker_type = linker_type

    def custom_reaction(self, reactants, rxn_smarts):
        '''
        This function performs a custom reaction on a set of reactants and SMARTS string.
        :param reactants:
        :param rxn_smarts:
        :return:
        '''
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        product = rxn.RunReactants(reactants)[0][0]
        Chem.SanitizeMol(product)
        return product

    def linker_to_ps(self, linker, ps) -> Chem.Mol:
        '''
        This function performs the reaction of the linker to the porphyrin.
        :param linker:
        :param ps:
        :return:
        '''
        if self.linker_type == 'NN':
            rxn_smarts = self.rcx_linker_to_ps_NN
        if self.linker_type == 'OO':
            rxn_smarts = self.rcx_linker_to_ps_OO
        return self.custom_reaction((linker, ps), rxn_smarts)

    def ligand_to_linked_ps(self, linked_ps, ligand) -> Chem.Mol:
        '''
        This function performs the reaction of the ligand to the linked porphyrin.
        :param linked_ps:
        :param ligand:
        :return:
        '''
        if self.linker_type == 'NN':
            rxn_smarts = self.rcx_linked_add_ligand_NN
        if self.linker_type == 'OO':
            rxn_smarts = self.rcx_linked_add_ligand_OO
        return self.custom_reaction((linked_ps, ligand), rxn_smarts)

    def conjugate(self, ps: Chem.Mol, linker: Chem.Mol, ligand: Chem.Mol) -> Chem.Mol:
        '''
        This function performs the entire conjugation reaction.
        :param ps:
        :param linker:
        :param ligand:
        :return:
        '''
        try:
            if type(ps) == str and type(linker) == str and type(ligand) == str:
                self.mode = 'smiles'
            if self.mode == 'smiles':
                ps, linker, ligand = self._fix_smiles(ps), self._fix_smiles(linker), self._fix_smiles(ligand)
                ps, linker, ligand = Chem.MolFromSmiles(ps), Chem.MolFromSmiles(linker), Chem.MolFromSmiles(ligand)
            self.set_linker_type(self._detect_linker_type(linker))
            linked_ps = self.linker_to_ps(linker, ps)
            conjugate = self.ligand_to_linked_ps(linked_ps, ligand)
            return conjugate
        except:
            # print('Error in conjugation reaction. Check input smiles.')
            # print('ps: ', Chem.MolToSmiles(ps))
            # print('linker: ', Chem.MolToSmiles(linker))
            # try:
            #     print('ligand: ', Chem.MolToSmiles(ligand))
            # except:
            #     print('ligand: ', ligand)
            return None

    def _detect_linker_type(self, linker: Chem.Mol) -> str:
        '''
        This function detects the type of linker used in the reaction.
        :param linker:
        :return:
        '''
        smiles = Chem.MolToSmiles(linker)
        if smiles[0] == 'N' and smiles[-1] == 'N':
            return 'NN'
        if smiles[0] == 'O' and smiles[-1] == 'O':
            return 'OO'

    def _fix_smiles(self, string):
        return repr(string)[1:-1]


if __name__ == '__main__':
    ps = r"OC(C1=CC=C(/C2=C3C=CC(/C(C4=C(C=CC=C4F)F)=C5N/C(C=C\5)=C(C6=C(C=CC=C6F)F)\C7=N/C(C=C7)=C(C8=C(C=CC=C8F)F)\C9=CC=C2N9)=N\3)C=C1)=O"
    linker_N = 'NCCCCCCN'
    linker = 'OCCCCCCO'
    ligand = 'Cc1c(c2cc(ccc2n1C(=O)c3ccc(cc3)Cl)OC)CC(=O)O'
    # reactants = [Chem.MolFromSmiles(x) for x in (ps, linker, ligand)]
    #
    # product = Reactor().conjugate(*reactants)
    # print(Chem.MolToSmiles(product))
    # Chem.MolToSmiles(product, kekuleSmiles=False)
    # #save the molecule as a sdf
    # Chem.MolToMolFile(product, 'product_test.sdf')
    # with open('/home/anton/PycharmProjects/cadd/data/parts/ps/ps.smi') as f:
    #     ps = f.read().splitlines()
    #     for i in ps:
    #         original_str = repr(i)[1:-1]
    #
    #         reactants = [Chem.MolFromSmiles(x) for x in (original_str, linker_N, ligand)]
    #         # print(F'Reading: {Chem.MolToSmiles(ps)}')
    #         product = Reactor().conjugate(*reactants)
    #         print(F'Writing:  {Chem.MolToSmiles(product)}')
    product = Reactor().conjugate(ps, linker_N, ligand)
    print(F'Writing:  {Chem.MolToSmiles(product)}')
