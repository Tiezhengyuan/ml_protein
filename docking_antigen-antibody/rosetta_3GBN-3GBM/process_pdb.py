'''
process PDB
1. renumber chains
'''
import os
from Bio.PDB import PDBIO, PDBParser

class ProcessPdb:

    def __init__(self, pdb_file:str):
        self.pdb_file = pdb_file
        self.pdb_file_name = os.path.basename(pdb_file)
        self.indir = os.path.dirname(pdb_file)

    def load_structure(self):
        """
        update self.structure
        """
        structure_id = self.pdb_file_name.replace('.pdb', '')
        with open(self.pdb_file, 'r') as f:
            parser = PDBParser(PERMISSIVE=1)
            self.structure = parser.get_structure(structure_id, f)

    def renumber(self, start:int=None, norestart=True, preserve=False):
        '''
        renumber chains
        arguments:
            norestart: don't start renumbering at each chain
            preserve: preserve insertion code and heteroflags
        '''
        residue_id = 1 if start is None else start

        # process structure
        chain_id = ""
        for residue in self.structure.get_residues():
            chain = residue.get_parent()
            if chain_id != chain.get_id() and not norestart:
                chain_id = chain.get_id()
                residue_id = int(start)

            if preserve:
                hetero = residue.id[0]
                insert = residue.id[2]
                residue.id = (hetero, residue_id, insert)
            else:
                residue.id = (' ', residue_id, ' ')
            residue_id += 1
            # print(chain.get_id(), end=',') 

    def get_table(self):
        """
        get the score table from the bottom of a PDB
        return a list of lines
        """
        raw_table = []
        tag = False
        with open(self.pdb_file, 'r') as f:
            for line in f:
                if line.startswith("#BEGIN_POSE_ENERGIES_TABLE"):
                    tag = True
                if len(line) > 1 and tag is True:
                    raw_table.append(line)
        return raw_table

    def export_pdb(self, outfile_name:str, outdir:str=None, has_table=False):
        outfile = os.path.join(outdir, outfile_name) if \
            outdir else os.path.join(self.indir, outfile_name)
        print(f"Try to create {outfile}...")
        with open(outfile, 'w') as f:
            # save atoms
            io=PDBIO()
            io.set_structure(self.structure)
            io.save(f)

            # save table
            if(has_table):
                raw_table = self.get_table()
                f.writelines(raw_table)