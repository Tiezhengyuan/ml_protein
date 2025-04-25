'''
build dataset from pdb
'''
import os
import numpy as np

alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
alpha_3 = [
    'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','GAP'
]

DEFAULT_ALPHABET = [chr(i) for i in range(65, 91)] + \
    [chr(i) for i in range(97, 123)] + \
    [str(item) for item in list(np.arange(300))]

class ParsePdb:

    def __init__(self, pdb_file:str, chain_list:list=None, ca_only:bool=False):
        self.pdb_file = pdb_file
        self.chain_list = chain_list if chain_list else DEFAULT_ALPHABET
        self.ca_only = ca_only
        self.pdb_dict = {
            'name': os.path.splitext(os.path.basename(self.pdb_file))[0],
        }

    def parse_PDB(self):
        '''
        parse PDB
        '''
        s = 0
        concat_seq = ''
        concat_N = []
        concat_CA = []
        concat_C = []
        concat_O = []
        concat_mask = []
        coords_dict = {}
        for letter in self.chain_list:
            xyz, seq = self.parse_PDB_biounits(letter)
            if type(xyz) != str:
                concat_seq += seq[0]
                self.pdb_dict['seq_chain_'+letter] = seq[0]
                coords_dict_chain = {}
                if self.ca_only:
                    coords_dict_chain['CA_chain_'+letter]=xyz.tolist()
                else:
                    coords_dict_chain['N_chain_' + letter] = xyz[:, 0, :].tolist()
                    coords_dict_chain['CA_chain_' + letter] = xyz[:, 1, :].tolist()
                    coords_dict_chain['C_chain_' + letter] = xyz[:, 2, :].tolist()
                    coords_dict_chain['O_chain_' + letter] = xyz[:, 3, :].tolist()
                self.pdb_dict['coords_chain_'+letter]=coords_dict_chain
                s += 1
        self.pdb_dict['num_of_chains'] = s
        self.pdb_dict['seq'] = concat_seq
        return self.pdb_dict


    def parse_PDB_biounits(self, chain=None):
        '''
        args:x = PDB filename
        args:atoms = atoms to extract (optional)
        return: (length, atoms, coords=(x,y,z)), sequence
        '''
        sidechain_atoms = ['CA'] if self.ca_only else ['N', 'CA', 'C', 'O']
        aa_3_N = {a:n for n,a in enumerate(alpha_3)}

        # check the best min_resn and max_resn
        xyz, seq, min_resn, max_resn = {}, {}, 1e6, -1e6
        for line in open(self.pdb_file,"rb"):
            line = line.decode("utf-8", "ignore").rstrip()

            if line[:6] == "HETATM" and line[17:17+3] == "MSE":
                line = line.replace("HETATM","ATOM  ")
                line = line.replace("MSE","MET")

            if line[:4] == "ATOM":
                ch = line[21:22]
                if ch == chain or chain is None:
                    atom = line[12:12+4].strip()
                    resi = line[17:17+3]
                    resn = line[22:22+5].strip()
                    x,y,z = [float(line[i:(i+8)]) for i in [30,38,46]]

                    if resn[-1].isalpha(): 
                        resa,resn = resn[-1],int(resn[:-1])-1
                    else: 
                        resa,resn = "",int(resn)-1
                    
                    # resn = int(resn)
                    if resn < min_resn: 
                        min_resn = resn
                    if resn > max_resn: 
                        max_resn = resn
                    if resn not in xyz: 
                        xyz[resn] = {}
                    if resa not in xyz[resn]: 
                        xyz[resn][resa] = {}
                    if resn not in seq: 
                        seq[resn] = {}
                    if resa not in seq[resn]: 
                        seq[resn][resa] = resi
                    if atom not in xyz[resn][resa]:
                        xyz[resn][resa][atom] = np.array([x,y,z])

        # convert to numpy arrays, fill in missing values
        seq_, xyz_ = [], []
        try:
            for resn in range(min_resn, max_resn+1):
                if resn in seq:
                    for k in sorted(seq[resn]):
                        seq_.append(aa_3_N.get(seq[resn][k],20))
                else:
                    seq_.append(20)
                if resn in xyz:
                    for k in sorted(xyz[resn]):
                        for atom in sidechain_atoms:
                            if atom in xyz[resn][k]:
                                xyz_.append(xyz[resn][k][atom])
                            else:
                                xyz_.append(np.full(3,np.nan))
                else:
                    for atom in sidechain_atoms:
                        xyz_.append(np.full(3, np.nan))
            xyz = np.array(xyz_).reshape(-1,len(sidechain_atoms),3)
            seq = self.N_to_AA(np.array(seq_))
            return xyz, seq 
        except TypeError:
            return 'no_chain', 'no_chain'

    def AA_to_N(self, x):
        states = len(alpha_1)
        aa_1_N = {a:n for n,a in enumerate(alpha_1)}
        # ["ARND"] -> [[0,1,2,3]]
        x = np.array(x);
        if x.ndim == 0:
            x = x[None]
        return [[aa_1_N.get(a, states-1) for a in y] for y in x]

    
    def N_to_AA(self, x):
        aa_N_1 = {n:a for n,a in enumerate(alpha_1)}
        # [[0,1,2,3]] -> ["ARND"]
        x = np.array(x);
        if x.ndim == 1:
            x = x[None]
        return ["".join([aa_N_1.get(a,"-") for a in y]) for y in x]