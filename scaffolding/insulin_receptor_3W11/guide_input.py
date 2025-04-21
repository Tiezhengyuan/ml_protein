'''
make guided inputs 
1. secondary structure *_ss.pt
2. block adjacency *_adj.pt
'''
from copy import deepcopy
import os
import random
import numpy as np
import pandas as pd
import torch
from Bio.PDB import PDBParser

from estimate_sse import EstimateSSE



class GuideInput:
    def __init__(self, atoms:list, outfile_prefix:str):
        '''
        args: atoms are list of residues of a certain chain
        '''
        self.atoms = atoms
        self.outfile_prefix = outfile_prefix
        # tensor of secondary structure
        self.ss = None
        # tensor of adjacency block
        self.adj = None

    def create_adj(self, cutoff:int=6, exclude_loops:bool=True):
        """
        convert adjacency to tensor dataset
        update self.adj
        """
        # initialize N14 representation
        n14_xyz = self.init_n14()
        # get X,Y,Z
        xyz = self.parse_xyz(n14_xyz)

        # estimate SSE
        sse = EstimateSSE(xyz[:,1]).get_sse()

        # convert SSE to int
        self.dummy_ss(sse)

        # 
        dist = self.get_dist(torch.tensor(xyz).float())

        # segments
        segments = self.get_segments(exclude_loops)

        # create block_adj tensor
        self.create_block_adj(dist, segments, cutoff)
        torch.save(self.adj, f'{self.outfile_prefix}_adj.pt')
        return self.adj

    def get_dist(self, xyz):
        '''
        estimate distance of atoms of N, alpha-C, and C
        NOTE: some hard coding
        '''
        MAX_DIST = 999.9

        # three anchor atoms
        N, Ca, C  = xyz[:,0], xyz[:,1], xyz[:,2]
        # generate_Cbeta
        # recreate Cb given N,Ca,C
        b = Ca - N 
        c = C - Ca
        # cross production along the last dim
        a = torch.cross(b, c, dim=-1)
        #Cb = -0.58273431*a + 0.56802827*b - 0.54067466*c + Ca
        # fd: below matches sidechain generator (=Rosetta params)
        Cb = -0.57910144 * a + 0.5689693 * b - 0.5441217 * c + Ca
    
        # May need a batch dimension - NRB
        dist = torch.cdist(Cb, Cb, p=2) # [L,L]
        dist[torch.isnan(dist)] = MAX_DIST
        dist += MAX_DIST * torch.eye(xyz.shape[0], device=xyz.device)
        return dist


    def get_segments(self, exclude_loops:bool) -> list:
        '''
        given self.ss, construct a list of segments 
        and the index at which they begin and end
        '''
        segments = []
        in_segment = True
        begin, end = -1, -1
        ss_num = len(self.ss)
        for i in range(ss_num):
            # Starting edge case
            if i == 0:
                begin = 0
            else:
                if not self.ss[i] == self.ss[i-1]:
                    end = i 
                    segments.append( (self.ss[i-1], begin, end) )
                    begin = i
        # Ending edge case: last segment is length one
        if not end == ss_num:
            segments.append( (self.ss[-1], begin, ss_num) )
        
        # filter segments
        if exclude_loops:
            segments = [s for s in segments if s[0] != 2]
        print('###segments;', segments)
        return segments

    def create_block_adj(self, dist, segments:list, cutoff:int):
        '''
        block adjacency atoms given a cutoff
        '''
        seg_num = len(segments)
        self.adj = torch.zeros_like(dist)
        for i in range(seg_num):
            curr_segment = segments[i]
            begin_i, end_i = curr_segment[1:3]
            for j in range(i+1, seg_num):
                j_segment = segments[j]
                begin_j, end_j = j_segment[1:3]
                if torch.any( dist[begin_i:end_i, begin_j:end_j] < cutoff ):
                    _ones = torch.ones(end_i - begin_i, end_j - begin_j)
                    # block_adj is symmetic
                    self.adj[begin_i:end_i, begin_j:end_j] = _ones
                    self.adj[begin_j:end_j, begin_i:end_i] = _ones

    def create_ss(self, min_mask:float=None, max_mask:float=None):
        """
        convert SSE to tensor dataset
        update self.ss
        """
        # initialize N14 representation
        n14_xyz = self.init_n14()
        # get X,Y,Z
        xyz = self.parse_xyz(n14_xyz)

        # estimate SSE
        sse = EstimateSSE(xyz[:,1]).get_sse()

        # convert SSE to int
        self.dummy_ss(sse)
        
        # to one-hot table
        self.ss = torch.nn.functional.one_hot(self.ss, num_classes=4)

        # update mask
        self.add_mask(min_mask, max_mask)

        # return ss_argmax for output
        ss_argmax = torch.argmax(self.ss, dim=1).float()
        torch.save(ss_argmax, f'{self.outfile_prefix}_ss.pt')
        return ss_argmax

    def init_n14(self, infile=None) -> dict:
        '''
        initialize N14 representation
        '''
        infile = infile if infile else 'AA_N14.csv'
        n14 = pd.read_csv(infile, header=0)
        n14 = n14.fillna('')
        n14 = n14.map(lambda a: a.strip())

        n14_xyz = {}
        for aa in list(n14):
            xyz = np.zeros((len(n14), 3), dtype='float')
            n14_xyz[aa] = pd.DataFrame(xyz, index=n14[aa], columns=list('xyz'))
        return n14_xyz

    def parse_xyz(self, n14_xyz):
        '''
        get X,Y,Z by atoms
        '''
        n14_len = np.max([len(i) for i in n14_xyz.values()])
        xyz = np.zeros((len(self.atoms), n14_len, 3))
        for i, res in enumerate(self.atoms):
            res_id = int(res.id[1])
            res_name = res.get_resname()
            # note: hard copy
            res_xyz = deepcopy(n14_xyz[res_name])
            _xyz = self.residue_xyz(res)
            # merge
            inter_names = list(set(res_xyz.index).intersection(_xyz.index))
            res_xyz.loc[inter_names] = _xyz.loc[inter_names]
            # update
            xyz[i] = res_xyz
        return xyz

    def residue_xyz(self, residue):
        '''
        get X,Y,Z of atoms given a residue
        '''
        xyz = {}
        for a in residue:
            xyz[a.id] = a.get_coord()
        xyz = pd.DataFrame(xyz, index=list('xyz')).T
        return xyz

    def dummy_ss(self, sse):
        '''
        convert character to integer
        '''
        ss_conv = {
            'H':0,  # helix
            'E':1,  # strand
            'L':2,  # loop
            'U':3,  #mask/unkown
        }
        self.ss = torch.tensor([int(ss_conv[i]) for i in sse])
        print('####ss:', self.ss)


    def add_mask(self, min_mask:float=None, max_mask:float=None):
        '''
        Add mask==3 into end of SSE block in self.ss 
        '''
        min_mask = 0 if min_mask is None else min_mask
        max_mask = 1.0 if max_mask is None else max_mask
        mask_prop = random.uniform(min_mask, max_mask)

        #gets last index of each block of ss
        transitions = np.where(self.ss[:-1] - self.ss[1:] != 0)[0]
        print('## transitions:', transitions)

        stuck_counter = 0
        unknown_percent = 0
        while unknown_percent < mask_prop or stuck_counter > 100:
            width = random.randint(1, 9)
            start = random.choice(transitions)
            offset = random.randint(-8, 1)
            start += offset
            end = start + width
            try:
                ss[start:end] = 3
            except:
                stuck_counter += 1
                pass
            unknown_percent = len(self.ss[self.ss == 3])/len(self.ss)

    def get_idx(self):
        '''
        '''
        idx = [int(res.id[1]) for res in self.atoms]
        return np.array(idx)

    def get_mask(self):
        '''
        get mask only from self.ss 
        NOTE: run after self.create_ss()
        '''
        idx = [int(res.id[1]) for res in self.atoms]

        # append residueNo to the last dimension
        mask_ss = torch.cat((self.ss, torch.tensor(idx)[..., None]), dim=-1)
        mask = np.argmax(mask_ss[:,:-1].numpy(), axis=-1)
        return torch.tensor(np.where(mask == 3))
