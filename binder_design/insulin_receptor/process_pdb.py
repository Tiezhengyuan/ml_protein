'''
process PDB
1. renumber chains
'''
import os
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import warnings

from amino_acids import longer_names


class NonHetSelect(Select):
    """Select only ATOM records (exclude HETATM)"""
    def accept_residue(self, residue):
        return residue.id[0] == " "  # Keep only standard residues

class ProcessPdb:

    def __init__(self, indir:str):
        self.indir = indir

    def load_structure(self, pdb_file_name:str):
        """
        update self.structure
        """
        self.pdb_file = os.path.join(self.indir, pdb_file_name)
        self.structure_id = pdb_file_name.replace('.pdb', '')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            parser = PDBParser(PERMISSIVE=1)
            self.structure = parser.get_structure(self.structure_id, self.pdb_file)
    
    def get_chains(self):
        for model in self.structure:
            for chain in model:
                print(f"Chain={chain.id}, AA residues={len(chain)}")

    def remove(self, chain_id:str, id_start:int, id_end:int):
        '''
        remove residues from a chain
        '''
        for model in self.structure:
            for chain in model:
                if chain.id == chain_id:
                    to_remove = [res for res in chain if id_start <= int(res.id[1]) <= id_end ]
                    # print(to_remove)
                    for res in to_remove:
                        chain.detach_child(res.id)

    def trim(self, chain_id:str, id_start:int):
        '''
        remove residues from a chain
        '''
        for model in self.structure:
            for chain in model:
                if chain.id == chain_id:
                    to_remove = [res for res in chain if int(res.id[1]) >= id_start ]
                    # print(to_remove)
                    for res in to_remove:
                        chain.detach_child(res.id)

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
    
    def export_fasta(self, outfile_name:str, outdir:str=None):
        '''
        retrieve sequence by chains from pdb
        '''
        outfile = os.path.join(outdir, outfile_name) if \
            outdir else os.path.join(self.indir, outfile_name)
        print(f"Try to create {outfile}...")

        records = []
        for model in self.structure:
            for chain in model:
                seq = []
                for residue in chain:
                    if residue.resname in longer_names:
                        seq.append(longer_names[residue.resname])
                rec = SeqRecord(
                    Seq(''.join(seq)),
                    id=f"{self.structure_id}|Chain {chain.id}",
                    description='',
                )
                records.append(rec)
        # save to fasta
        with open(outfile, 'w') as f:
            SeqIO.write(records, f, 'fasta')

    def split_by_chain(self, chain_id:str):
        """
            Save chain while preserving original header information
        """
        # Copy header from original file
        header_lines = []
        with open(self.pdb_file) as f:
            for line in f:
                if line.startswith(("HEADER", "TITLE", "COMPND", "SOURCE")):
                    header_lines.append(line)

        # Write output with header
        io = PDBIO()
        io.set_structure(self.structure[0][chain_id])

        outfile = os.path.join(self.indir, f"{self.structure_id}_{chain_id}.pdb")        
        with open(outfile, "w") as f:
            f.writelines(header_lines)
            io.save(f, NonHetSelect())
        print(f"Saved chain {chain_id} with metadata to {outfile}")


