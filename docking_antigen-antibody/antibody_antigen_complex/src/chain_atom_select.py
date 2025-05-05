'''
'''
from Bio.PDB import Select, PDBIO

class ChainAtomSelect(Select):
    def __init__(self, chains):
        self.chains = chains
        
    def accept_chain(self, chain):
        if chain.get_id() in self.chains:
            return True
        return False
        
    def accept_residue(self, residue):
        if residue.get_id()[0] == ' ':
            return True
        return False