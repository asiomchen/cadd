from rdkit import Chem
from reaction import Reactor
import unittest


class TestReactor(unittest.TestCase):
    def setUp(self):
        self.reactor = Reactor()
        self.acid = 'CCC(O)=O'
        self.amine = 'NCCCCN'
        self.acid = Chem.MolFromSmiles(self.acid)
        self.amine = Chem.MolFromSmiles(self.amine)
        self.linked = self.reactor.linker_to_ps(self.amine, self.acid)
        self.conjugate = self.reactor.ligand_to_linked_ps(self.linked, self.acid)
        self.whole = self.reactor.conjugate(self.acid, self.amine, self.acid)

    def test_linker_to_ps(self):
        self.assertEqual(Chem.MolToSmiles(self.linked), 'CCC(=O)NCCCCN')

    def test_ligand_to_linked_ps(self):
        self.assertEqual(Chem.MolToSmiles(self.conjugate), 'CCC(=O)NCCCCNC(=O)CC')

    def test_conjunction(self):
        self.assertEqual(Chem.MolToSmiles(self.whole), 'CCC(=O)NCCCCNC(=O)CC')
