'''
Split a protein-ligand complex into protein and ligands.
Assign ligand bond orders using SMILES strings from Ligand Export
'''

from functools import partial
import gzip
import re
import sys
from io import StringIO
import requests
import numpy as np
import pandas as pd
import os
import random
import traceback
import torch
from pprint import pprint

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import prody
import openfold
from openfold.data import data_transforms
from openfold.np import protein, residue_constants

from constants import punctuation_regex, molecule_regex, ubiquitous_ligands


class ParseComplex:

    def __init__(self, pdb_path_file, mol_wt_cutoff=None, min_atoms=None):
        self.pdb_path_file = pdb_path_file
        self.data_dir = os.path.dirname(pdb_path_file)
        # minimum molecular weight of a ligand >= 100Da
        self.mol_wt_cutoff = mol_wt_cutoff if mol_wt_cutoff else 100
        # minimum atoms of a ligand >= 3 residues
        self.min_atoms = min_atoms if min_atoms else 3
        # 
        self.ligand_expo = None

    def parse(self, nrows:int=None):
        name = os.path.splitext(os.path.basename(self.pdb_path_file))[0]
        with open(self.pdb_path_file, 'r') as f:
            files = f.read().split('\n')
        if nrows:
            files = files[:nrows]

        # chunksize dataset
        pool = []
        for pdb_file in files:
            res = self.process_entry(pdb_file)
            if res:
                pool.extend(res)
        print(len(pool), ' pdb files are parsed.')
        
        # export
        if pool:
            df = pd.DataFrame(pool)
            outfile = os.path.join(self.data_dir, f'{name}.p')
            df.to_pickle(outfile)

            # delete 
            os.remove(self.pdb_path_file)



################################################33
    # def split_pdb_files(self, train_size=.9):
    #     pdb_files = [i for i in self.iter_pdb_file()]
    #     pdb_num = len(pdb_files)
    #     print(f"total pdb files: {pdb_num}")
        
    #     # shuffle
    #     pdb_files = sorted(pdb_files)
    #     random.seed(42)
    #     random.shuffle(pdb_files)
        
    #     # split
    #     split_idx = int(train_size * pdb_num)
    #     train = pdb_files[:split_idx]
    #     test = pdb_files[split_idx:]
    #     return train, test
    
    # def iter_pdb_file(self):
    #     for root, dirs, names in os.walk(self.pdb_dir):
    #         for file_name in names:
    #             if file_name.endswith('.gz'):
    #                 yield os.path.join(root, file_name)

    def read_ligand_expo(self):
        """
        Get Ligand Expo data from http://ligand-expo.rcsb.org/
        Ligand Expo (formerly Ligand Depot) provides chemical 
        and structural information about small molecules 
        within the structure entries of the Protein Data Bank. 
        """
        # try to find a file Components-smiles-stereo-oe.smi
        file_name = "Components-smiles-stereo-oe.smi"
        infile = os.path.join(self.data_dir, file_name)
        if not os.path.isfile(infile):
            url = f"http://ligand-expo.rcsb.org/dictionaries/{file_name}"
            print(f"Try to download data from RCSB: {url}, and save them into {infile}.")
            r = requests.get(url, allow_redirects=True)
            open(infile, 'wb').write(r.content)

        # load data to self.ligand_expo
        self.ligand_expo = pd.read_csv(
            infile,
            sep="\t",
            header=None,
            names=["SMILES", "ID", "Name"],
            na_filter=False
        )


    def process_entry(self, pdb_file):
        """
        Slit pdb into protein and ligands,
        parse protein sequence and ligand tokens
        """
        protein, ligand = self.get_prody_structure(pdb_file)
        if ligand is None or protein is None:
            return []

        # get protein
        chains = [chain.getChid() for chain in protein.getHierView()]
        pdb_str = self.get_pdb_str(protein)
        seq, features = self.get_protein_sequence_and_coords(chains, pdb_str)
        if seq is None or features is None:
            return []

        entries = []
        # filter ligands by molecular weight
        res_names = []
        for i in ligand.getResnames():
            if i not in ubiquitous_ligands and i not in res_names:
                res_names.append(i)
        for res_name in res_names:
            mols, template = self.process_ligand(ligand, res_name)
            if mols and template:
                mol_wt = Descriptors.ExactMolWt(template)
                natoms = template.GetNumAtoms()
                if mol_wt >= self.mol_wt_cutoff and natoms >= self.min_atoms:
                    # only use first copy of ligand
                    ligand_mol = mols[0]
                    # smiles and xyz
                    smi, xyz, xyz_2d, atom_idx = self.tokenize_ligand(ligand_mol)
                    # ligand bonds
                    bonds, ligand_token_bonds = [], []
                    for b in template.GetBonds():
                        bond = (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
                        bonds.append(bond)
                        token_bond = (atom_idx.index(bond[0]), atom_idx.index(bond[1]))
                        ligand_token_bonds.append(token_bond)
                    # update
                    entry = {
                        'pdb_id': os.path.basename(pdb_file).split('.')[-3][3:],
                        'seq': seq,
                        'rigidgroups_gt_frames': features['rigidgroups_gt_frames'],
                        'torsion_angles_sin_cos': features['torsion_angles_sin_cos'],
                        'lig_id': res_name,
                        'ligand_mol': ligand_mol,
                        'smiles': smi,
                        'ligand_xyz': xyz,
                        'ligand_xyz_2d': xyz_2d,
                        'ligand_bonds': bonds,
                        'ligand_token_bonds': ligand_token_bonds,
                    }
                    entries.append(entry)
        print(len(entries))
        return entries

    def get_prody_structure(self, pdb_file):
        """
        Split a protein-ligand pdb into protein and ligand components
        """
        with gzip.open(pdb_file,'rt') as f:
            pdb = prody.parsePDBStream(f)
        print(pdb_file)
        protein = pdb.select('protein')
        ligand = pdb.select('not protein and not water')
        return protein, ligand

    def process_ligand(self, ligand, res_name):
        """
        Add bond orders to a pdb ligand
        1. Select the ligand component with name "res_name"
        2. Get the corresponding SMILES from the Ligand Expo dictionary
        3. Create a template molecule from the SMILES in step 2
        4. Write the PDB file to a stream
        5. Read the stream into an RDKit molecule
        6. Assign the bond orders from the template from step 3
        :param ligand: ligand as generated by prody
        :param res_name: residue name of ligand to extract
        :return: molecule with bond orders assigned
        """
        df = self.ligand_expo
        sub_smiles = df[df['ID'].values == res_name]['SMILES'].values[0]
        template = AllChem.MolFromSmiles(sub_smiles)

        allres = ligand.select(f"resname {res_name}")
        res = np.unique(allres.getResindices())
        mols = []
        for i in res:
            # print(f"resname {res_name} and resindex {i}")
            sub_mol = ligand.select(f"resname {res_name} and resindex {i}")
            output = StringIO()
            prody.writePDBStream(output, sub_mol)
            pdb_string = output.getvalue()
            rd_mol = AllChem.MolFromPDBBlock(pdb_string)
            try:
                mol = AllChem.AssignBondOrdersFromTemplate(template, rd_mol)
                mols.append(mol)
            except Exception as e:
                # print(f"ligand={res_name}, error={e}")
                pass
        return mols, template

    def get_pdb_str(self, protein):
        # structure of openfold
        pdb_structure = StringIO()
        prody.writePDBStream(pdb_structure, protein)
        pdb_str = pdb_structure.getvalue()
        return pdb_str

    def get_protein_sequence_and_coords(self, chains, pdb_str):
        '''
        chains:
        pdb_str: string type, atoms lines in pdb
        '''
        aatype = []
        atom_positions = []
        atom_mask = []
        for chain in chains:
            try:
                _p = None
                _p = protein.from_pdb_string(pdb_str, chain)
                aatype.append(_p.aatype)
                atom_positions.append(_p.atom_positions)
                atom_mask.append(_p.atom_mask)
            except Exception as e:
                print(pdb_str[:10], chain, e)
        
        if len(aatype) == len(atom_positions) == len(atom_mask) >0:
            # concatenate chains
            aatype = np.concatenate(aatype)
            atom_positions = np.concatenate(atom_positions)
            atom_mask = np.concatenate(atom_mask)
            
            # sequence
            seq = residue_constants.aatype_to_str_sequence(aatype)

            # determine torsion angles
            features = {
                'aatype': torch.tensor(aatype),
                'all_atom_positions': torch.tensor(atom_positions),
                'all_atom_mask': torch.tensor(atom_mask)
            }
            features = data_transforms.atom37_to_torsion_angles()(features)
            features = data_transforms.atom37_to_frames(features)
            features = data_transforms.make_atom14_masks(features)
            features = data_transforms.make_atom14_positions(features)
            features = {k: v.numpy() for k, v in features.items() \
                if isinstance(v, torch.Tensor)}
            return seq, features
        return None, None

    def tokenize_ligand(self, mol):
        # convert to SMILES and map atoms
        smi = Chem.MolToSmiles(mol)

        # position of atoms in SMILES (not counting punctuation)
        atom_order = [int(s) for s in list(filter(None,re.sub(r'[\[\]]','',mol.GetProp("_smilesAtomOutputOrder")).split(',')))]

        # tokenize the SMILES
        tokens = list(filter(None, re.split(molecule_regex, smi)))

        # remove punctuation
        masked_tokens = [re.sub(punctuation_regex,'',s) for s in tokens]

        k = 0
        token_pos = []
        for i,token in enumerate(masked_tokens):
            if token != '' and k < len(atom_order):
                pos = mol.GetConformer().GetAtomPosition(atom_order[k])
                token_pos.append(tuple(pos))
                k += 1
            else:
                token_pos.append((np.nan, np.nan, np.nan))

        k = 0
        conf_2d = Chem.AllChem.Compute2DCoords(mol)
        token_pos_2d = []
        atom_idx = []
        for i,token in enumerate(masked_tokens):
            if token != '' and k < len(atom_order):
                pos_2d = mol.GetConformer(conf_2d).GetAtomPosition(atom_order[k])
                token_pos_2d.append(tuple(pos_2d))
                atom_idx.append(atom_order[k])
                k += 1
            else:
                token_pos_2d.append((0.,0.,0.))
                atom_idx.append(None)

        return smi, token_pos, token_pos_2d, atom_idx

if __name__ == '__main__':
    pdb_path_file = sys.argv[1]
    parser = ParseComplex(pdb_path_file)
    # read ligand table, update self.ligand_expo
    parser.read_ligand_expo()

    # 
    parser.parse()