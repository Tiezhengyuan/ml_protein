'''
process PDB
1. renumber chains
'''
import gzip
import os
from Bio.PDB import PDBIO, PDBParser, MMCIFParser, Select
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Filter out the specific warning
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter("ignore", PDBConstructionWarning)

from dir import Dir
from amino_acids import longer_names
from chain_atom_select import ChainAtomSelect


class NonHetSelect(Select):
    """Select only ATOM records (exclude HETATM)"""
    def accept_residue(self, residue):
        return residue.id[0] == " "  # Keep only standard residues

class ProcessPdb:

    def __init__(self, infile:str, outdir:str=None):
        # pdb_file
        self.pdb_file = infile
        self.indir = os.path.dirname(infile)
        self.pdb_file_name = os.path.basename(infile)
        self.outdir = outdir if outdir else self.indir
        # default objects should be updated
        self.chains = []
        # load structure
        self.get_parser()
        self.load_structure()

    def get_parser(self):
        '''
        update self.parser
        '''
        file_type = 'cif' if self.pdb_file_name.endswith('.cif') or \
            self.pdb_file_name.endswith('.cif.gzs')else 'pdb'

        self.parse = None
        if file_type == 'pdb':
            self.parser = PDBParser(PERMISSIVE=1)
        elif file_type == 'cif':
            self.parser = MMCIFParser(QUIET=True)

    def load_structure(self):
        """
        update self.structure
        """
        self.structure = None
        if self.parser:
            if self.pdb_file.endswith('.gz'):
                with gzip.open(self.pdb_file, 'rt') as gzf:
                    self.structure = self.parser.get_structure(self.pdb_file_name, gzf)
            else:
                self.structure = self.parser.get_structure(self.pdb_file_name, self.pdb_file)
            if self.structure:
                self.structure_id = self.structure.header.get('idcode', '').upper()
                print(f"Successfully retrieve structure of {self.structure_id}")
            # update self.outdir
            prefix = self.structure_id[:2]
            self.outdir = os.path.join(self.outdir, prefix, self.structure_id)
            Dir(self.outdir).init_dir()
            print('outputs dir: ', self.outdir)

    def get_header(self):
        header = {
            'resolution': self.structure.header['resolution'],
            'keywords': self.structure.header['keywords'],
            'compound': self.structure.header['compound'],
            'source': self.structure.header['source'],
        }
        return header

    def iter_res(self):
        '''
        iterate residues
        '''
        for model in self.structure:
            for chain in model:
                for res in chain:
                    yield model, chain, res

    def get_chains(self):
        for model in self.structure:
            _chains = {
                'model_id': model.id,
                'chains': [],
            }
            for chain in model:
                _chain = {
                    'chain_id': chain.id,
                    'AA residues': len(chain),
                    'seq': self.get_chain_seq(chain)
                }
                _chains['chains'].append(_chain)
            self.chains.append(_chains)
        return self.chains

    def get_chain_seq(self, chain):
        seq = []
        for residue in chain:
            if residue.resname in longer_names:
                seq.append(longer_names[residue.resname])
        return ''.join(seq)

    def parse_alpha_atoms(self, chain_id:str):
        '''
        Get all alpha atoms from Chain A
        '''
        for model in self.structure:
            for chain in model:
                if chain.id == chain_id:
                    ca_atoms = []
                    for residue in chain:
                        # Check if residue has a alpha atom
                        if "CA" in residue: 
                            ca_atoms.append(residue)
                    return ca_atoms
        return None

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

    def split_by_chain(self, chain_ids:list):
        """
            Save chain while preserving original header information
        """
        # Copy header from original file
        header_lines = []
        if self.pdb_file.endswith('.gz'):
            with gzip.open(self.pdb_file, 'rt') as f:
                for line in f:
                    if line.startswith(("HEADER", "TITLE", "COMPND", "SOURCE")):
                        header_lines.append(line)
        else:
            with open(self.pdb_file, 'r') as f:
                for line in f:
                    if line.startswith(("HEADER", "TITLE", "COMPND", "SOURCE")):
                        header_lines.append(line)

        # Write output with header
        io = PDBIO()
        io.set_structure(self.structure)
        
        outname = '_'.join([self.structure_id] + chain_ids) + '.pdb'
        outfile = os.path.join(self.outdir, outname)
        with open(outfile, "w") as f:
            f.writelines(header_lines)
            io.save(f, ChainAtomSelect(chain_ids))
        print(f"Saved chain {chain_ids} to {outfile}")
        return outfile


