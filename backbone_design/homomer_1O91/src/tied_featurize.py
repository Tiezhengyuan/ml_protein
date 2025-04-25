from __future__ import print_function
import json, time, os, sys, glob
import shutil
import numpy as np
import torch
from torch import optim
from torch.utils.data import DataLoader
from torch.utils.data.dataset import random_split, Subset

from copy import deepcopy
import torch.nn as nn
import torch.nn.functional as F
import random
import itertools
from Bio.Align import PairwiseAligner

from protein_mpnn_utils import _scores, _S_to_seq

class TiedFeaturize:

    def __init__(self, conf, protein):
        '''
        batch, device,and chain_dict should be declared
        '''
        self.conf = conf
        self.protein = protein

        self.device = conf.device_name
        # initialize
        self._chain_id_list_list = []
        self._visible_list_list = []
        self._masked_list_list = []
        self._masked_chain_length_list_list = []
        self._tied_pos_list_of_lists_list = []
        self._tied_beta = []


    def run(self):
        """
        Pack and pad batch into torch tensors
        """
        self.batch = [deepcopy(self.protein) for i in range(self.conf.batch_size)]
        # number of sequences in one batch
        B = self.conf.batch_size
        # max length of sequences in batch
        L_max = max([len(b['seq']) for b in self.batch])

        # X and S
        X_shape = [B, L_max, 1, 3] if self.conf.ca_only else [B, L_max, 4, 3]
        self._X = np.zeros(X_shape)
        self._S = np.zeros([B, L_max], dtype=np.int32)

        # other attributes
        self._residue_idx = -100 * np.ones([B, L_max], dtype=np.int32)
        #1.0 for the bits that need to be predicted
        self._chain_M = np.zeros([B, L_max], dtype=np.int32)
        #1.0 for the bits that need to be predicted
        self.pssm_coef_all = np.zeros([B, L_max], dtype=np.float32)
        #1.0 for the bits that need to be predicted
        self.pssm_bias_all = np.zeros([B, L_max, 21], dtype=np.float32)
        #1.0 for the bits that need to be predicted
        self._pssm_log_odds_all = 10000.0 * np.ones([B, L_max, 21], dtype=np.float32)
        #1.0 for the bits that need to be predicted
        self._chain_M_pos = np.zeros([B, L_max], dtype=np.int32)
        self._bias_by_res_all = np.zeros([B, L_max, 21], dtype=np.float32)
        #1.0 for the bits that need to be predicted
        self._chain_encoding_all = np.zeros([B, L_max], dtype=np.int32) 
        self._omit_AA_mask = np.zeros([B, L_max, len(self.conf.alphabet)], dtype=np.int32)
       
        # Build the batch
        self.initialize_chains()
        
        for i, batch in enumerate(self.batch):
            a, c = 0, 1
            l0, l1 = 0, 0
            mask_dict = {}
            visible_list = []
            masked_list = []
            masked_chain_length_list = []
            list_pool = {
                'global_idx_start': [0],
                'chain_seq': [],
                'x_chain': [],
                'chain_mask': [],
                'chain_encoding': [],
                'fixed_position_mask': [],
                'omit_AA_mask': [],
                'pssm_coef': [],
                'pssm_bias': [],
                'pssm_log_odds': [],
                'bias_by_res': [],
            }
            for step, chain_id in enumerate(self.all_chains):
                if chain_id in self.visible_chains:
                    visible_list.append(chain_id)
                    chain_seq, chain_length = self.update_chain_seq(list_pool, batch, chain_id)

                    # x_chain
                    self.update_x_chain(list_pool, batch, chain_id)
                    # chain mask
                    self.update_chain(list_pool, chain_id, chain_length, c)

                    l1 += chain_length
                    self._residue_idx[i, l0:l1] = 100 * (c-1) + np.arange(l0, l1)
                    l0 += chain_length
                    c += 1

                    # position mask
                    self.update_position_mask(list_pool, batch, chain_id, chain_length)
                    # AA mask
                    self.update_AA_mask(list_pool, batch, chain_id, chain_length)
                    # pssm
                    self.update_pssm(list_pool, batch, chain_id, chain_length)
                    # bias by res
                    self.update_bias_by_res(list_pool, batch, chain_id, chain_length)
                
                if chain_id in self.masked_chains:
                    masked_list.append(chain_id)
                    chain_seq, chain_length = self.update_chain_seq(list_pool, batch, chain_id)
                    masked_chain_length_list.append(chain_length)

                    # x_chain
                    self.update_x_chain(list_pool, batch, chain_id)
                    # chain mask
                    self.update_chain(list_pool, chain_id, chain_length, c)

                    l1 += chain_length
                    self._residue_idx[i, l0:l1] = 100 * (c-1) + np.arange(l0, l1)
                    l0 += chain_length
                    c += 1

                    # position mask
                    self.update_position_mask(list_pool, batch, chain_id, chain_length)
                    # AA_mask
                    self.update_AA_mask(list_pool, batch, chain_id, chain_length)
                    # pssm
                    self.update_pssm(list_pool, batch, chain_id, chain_length)
                    # bias by res
                    self.update_bias_by_res(list_pool, batch, chain_id, chain_length)
            
            # update self._tied_pos_list_of_lists_list
            # update self._tied_beta
            self.update_tied(L_max, batch, list_pool['global_idx_start'])
   
            # pad attributes then update global objects
            all_sequence = "".join(list_pool['chain_seq'])
            l = len(all_sequence)
            # x_chain: [L, 4, 3]
            self._X[i,:,:,:] = self.pad_attr(list_pool['x_chain'], [[0, L_max-l], [0, 0], [0, 0]])
            #[L,], 1.0 for places that need to be predicted
            self._chain_encoding_all[i,:] = self.pad_attr(list_pool['chain_encoding'], [[0, L_max-l]])
            self._chain_M[i,:] = self.pad_attr(list_pool['chain_mask'], [[0, L_max-l]])
            self._chain_M_pos[i,:] = self.pad_attr(list_pool['fixed_position_mask'], [[0, L_max-l]])
            self._omit_AA_mask[i,] = self.pad_attr(list_pool['omit_AA_mask'], [[0, L_max-l]])
            self.pssm_coef_all[i,:] = self.pad_attr(list_pool['pssm_coef'], [[0, L_max-l]])
            self.pssm_bias_all[i,:] = self.pad_attr(list_pool['pssm_bias'], [[0, L_max-l], [0, 0]])
            self._pssm_log_odds_all[i,:] = self.pad_attr(list_pool['pssm_log_odds'], [[0, L_max-l], [0, 0]])
            #[L,21], 0.0 for places where AA frequencies don't need to be tweaked
            self._bias_by_res_all[i,:] = self.pad_attr(list_pool['bias_by_res'], [[0, L_max-l], [0, 0]])

            # Convert to labels
            indices = np.asarray([self.conf.alphabet.index(a) for a in all_sequence], dtype=np.int32)
            self._S[i, :l] = indices
            self._visible_list_list.append(visible_list)
            self._masked_list_list.append(masked_list)
            self._masked_chain_length_list_list.append(masked_chain_length_list)

        # dihedral mask
        jumps = (
                (self._residue_idx[:,1:] - self._residue_idx[:,:-1]) == 1
            ).astype(np.float32)
        # shape [B,L,3]
        phi_mask = np.pad(jumps, [[0,0],[1,0]])[:,:,None]
        psi_mask = np.pad(jumps, [[0,0],[0,1]])[:,:,None]
        omega_mask = np.pad(jumps, [[0,0],[0,1]])[:,:,None]
        self._dihedral_mask = np.concatenate([phi_mask, psi_mask, omega_mask], axis=-1) 

        return None

    def initialize_chains(self):
        '''
        initialize chains
        '''
        self.all_chains, self.masked_chains, self.visible_chains = [], [], []
        for i, batch in enumerate(self.batch):
            if self.conf.chain_dict != None:
                #masked_chains a list of chain chain_ids to predict [A, D, F]
                self.masked_chains, self.visible_chains = self.conf.chain_dict[batch['name']]
            else:
                self.masked_chains = [item[-1:] for item in list(batch) if item[:10] == 'seq_chain_']
                self.visible_chains = []
            self.masked_chains.sort() #sort masked_chains 
            self.visible_chains.sort() #sort visible_chains 
            self.all_chains = self.masked_chains + self.visible_chains
        print(self.all_chains, self.masked_chains, self.visible_chains)


    def update_chain_seq(self, list_pool, batch, chain_id):
        '''
        update chain_seq
        '''
        chain_seq = batch[f'seq_chain_{chain_id}']
        chain_seq = ''.join([a if a!='-' else 'X' for a in chain_seq])
        chain_length = len(chain_seq)
        _val = list_pool['global_idx_start'][-1] + chain_length
        # update
        list_pool['global_idx_start'].append(_val)
        list_pool['chain_seq'].append(chain_seq)
        return chain_seq, chain_length

    def update_x_chain(self, list_pool, batch, chain_id):
        '''
        update x_chain
        '''
        chain_coords = batch[f'coords_chain_{chain_id}'] #this is a dictionary
        if self.conf.ca_only:
            #[chain_lenght,1,3] #CA_diff
            x_chain = np.array(chain_coords[f'CA_chain_{chain_id}'])
            if len(x_chain.shape) == 2:
                x_chain = x_chain[:,None,:]
            return x_chain

        chain_names = [f"{i}_chain_{chain_id}" for i in ['N', 'CA', 'C', 'O']]
        #[chain_lenght,4,3]
        x_chain = np.stack([chain_coords[c] for c in chain_names], 1) 
        # update _list
        list_pool['x_chain'].append(x_chain)


    def update_chain(self, list_pool, chain_id:str, chain_length:int, c:int):
        '''
        update chain_mask, chain_seq
        '''
        # chain mask
        #0.0 for visible chains, #1.0 for masked
        chain_mask = np.ones(chain_length) if chain_id in self.masked_chains \
            else np.zeros(chain_length) 
        list_pool['chain_mask'].append(chain_mask)

        # chain encoding
        chain_encoding = c * np.ones(np.array(chain_mask).shape[0])
        list_pool['chain_encoding'].append(chain_encoding)

    def update_position_mask(self, list_pool, batch, chain_id, chain_length):
        '''
        update fixed_position_mask
        '''
        fixed_position_mask = np.ones(chain_length)
        if chain_id in self.masked_chains and self.conf.fixed_position_dict:
                fixed_pos_list = self.conf.fixed_position_dict[batch['name']][chain_id]
                if fixed_pos_list:
                    fixed_position_mask[np.array(fixed_pos_list)-1] = 0.0
        list_pool['fixed_position_mask'].append(fixed_position_mask)

    def update_AA_mask(self, list_pool, batch, chain_id, chain_length):
        '''
        update omit_AA_mask
        '''
        omit_AA_mask_temp = np.zeros([chain_length, len(self.conf.alphabet)], np.int32)
        if chain_id in self.masked_chains and self.conf.omit_AA_dict:
            for item in self.conf.omit_AA_dict[batch['name']][chain_id]:
                idx_AA = np.array(item[0])-1
                AA_idx = np.array(
                    [np.argwhere(np.array(list(self.conf.alphabet))== AA)[0][0] for AA in item[1]]).repeat(idx_AA.shape[0]
                )
                idx_ = np.array([[a, b] for a in idx_AA for b in AA_idx])
                omit_AA_mask_temp[idx_[:,0], idx_[:,1]] = 1
        list_pool['omit_AA_mask'].append(omit_AA_mask_temp)

    def update_pssm(self, list_pool, batch, chain_id, chain_length):
        '''
        update pssm_bias, pssm_coef, pssm_log_odds
        '''
        # 
        pssm_coef = np.zeros(chain_length)
        pssm_bias = np.zeros([chain_length, 21])
        pssm_log_odds = 10000.0 * np.ones([chain_length, 21])
        if self.conf.pssm_dict:
            if chain_id in self.masked_chains and self.conf.pssm_dict[batch['name']][chain_id]:
                pssm_coef = self.conf.pssm_dict[batch['name']][chain_id]['pssm_coef']
                pssm_bias = self.conf.pssm_dict[batch['name']][chain_id]['pssm_bias']
                pssm_log_odds = self.conf.pssm_dict[batch['name']][chain_id]['pssm_log_odds']
        list_pool['pssm_coef'].append(pssm_coef)
        list_pool['pssm_bias'].append(pssm_bias)
        list_pool['pssm_log_odds'].append(pssm_log_odds)

    def update_bias_by_res(self, list_pool, batch, chain_id, chain_length):
        '''
        update bias_by_res_list
        '''
        bias_by_res = np.zeros([chain_length, 21])
        if chain_id in self.masked_chains and self.conf.bias_by_res_dict:
            bias_by_res = self.conf.bias_by_res_dict[batch['name']][chain_id]
        # update
        list_pool['bias_by_res'].append(bias_by_res)

    def pad_attr(self, _list:list, _dim:list):
        '''
        transform list object in list_pool
        '''
        # statck firstly
        _stack = np.concatenate(_list, 0)
        #default [L,], 1.0 for places that need to be predicted
        return np.pad(_stack, _dim, 'constant', constant_values=(0.0,))


    def update_tied(self, L_max, batch, global_idx_start_list):
        '''
        update self._chain_id_list_list
        update self._tied_pos_list_of_lists_list
        update self._tied_beta
        '''
        chain_id_list = []
        for chain_id in self.all_chains:
            if chain_id in self.visible_chains:
                chain_id_list.append(chain_id)
            if chain_id in self.masked_chains:
                chain_id_list.append(chain_id)

        tied_beta = np.ones(L_max)
        tied_pos_list_of_lists = []
        if self.conf.tied_positions_dict != None:
            chain_id_list_np = np.array(chain_id_list)
            tied_pos_list = self.conf.tied_positions_dict[batch['name']]
            if tied_pos_list:
                set_chains_tied = set(list(itertools.chain(*[list(item) for item in tied_pos_list])))
                for tied_item in tied_pos_list:
                    one_list = []
                    for k, v in tied_item.items():
                        start_idx = global_idx_start_list[np.argwhere(chain_id_list_np == k)[0][0]]
                        if isinstance(v[0], list):
                            for v_count in range(len(v[0])):
                                one_list.append(start_idx+v[0][v_count]-1)#make 0 to be the first
                                tied_beta[start_idx+v[0][v_count]-1] = v[1][v_count]
                        else:
                            for v_ in v:
                                one_list.append(start_idx+v_-1)#make 0 to be the first
                    tied_pos_list_of_lists.append(one_list)
        # update
        self._chain_id_list_list.append(chain_id_list)
        self._tied_pos_list_of_lists_list.append(tied_pos_list_of_lists)
        self._tied_beta = torch.from_numpy(tied_beta).to(dtype=torch.float32, device=self.device)




    def native_seq_info(self, log_probs):
        '''
        origin template seq
        '''
        native_seq = _S_to_seq(self.S[0], self.chain_M[0])
        native_score = _scores(self.S, log_probs, self.mask_for_loss)

        sorted_masked_chain_chain_ids = np.argsort(self.masked_list_list[0])
        print_masked_chains = [self.masked_list_list[0][i] for i in sorted_masked_chain_chain_ids]
        sorted_visible_chain_chain_ids = np.argsort(self.visible_list_list[0])
        print_visible_chains = [self.visible_list_list[0][i] for i in sorted_visible_chain_chain_ids]
        native_score_print = np.format_float_positional(
            np.float32(native_score.mean()), unique=False, precision=4
        )
        
        info = {
            'seq_type': 'native',
            'model_name': self.model_name,
            'name': self.name,
            'seq': native_seq,
            'masked_seq': self.process_mask_seq(native_seq, 0),
            'score': native_score,
            'score_print': native_score_print,
            'visible_chains': print_visible_chains,
            'designed_chains': print_masked_chains,
        }
        # print(info)
        return info

    def get_sample_seq(self, sample_dict) -> list:
        '''
        '''
        S_sample = sample_dict['S']
        res = []
        for b_ix in range(self.conf.batch_size):
            seq = _S_to_seq(S_sample[b_ix], self.chain_M[b_ix])
            seq = self.process_mask_seq(seq, b_ix)
            res.append(seq)
        return res


    def process_mask_seq(self, seq, b_ix):
        '''
        '''
        masked_chain_length_list = self.masked_chain_length_list_list[b_ix]
        masked_list = self.masked_list_list[b_ix]

        start = 0
        end = 0
        list_of_AAs = []
        for mask_l in masked_chain_length_list:
            end += mask_l
            list_of_AAs.append(seq[start:end])
            start = end
        seq = "".join(list(np.array(list_of_AAs)[np.argsort(masked_list)]))

        # mask
        # insert slash per mask length into seq
        l0 = 0
        _mask_chain = np.array(masked_chain_length_list)[np.argsort(masked_list)]
        for mc_length in list(_mask_chain)[:-1]:
            l0 += mc_length
            seq = seq[:l0] + '/' + seq[l0:]
            l0 += 1        
        return seq

    def cal_seq_recovery_rate(self, sample_dict) -> list:
        '''
        recovery rate: 0-1
        '''
        S_sample = sample_dict["S"]
        res = []
        for b_ix in range(self.conf.batch_size):
            _a = F.one_hot(self.S[b_ix], 21) * F.one_hot(S_sample[b_ix], 21)
            _b = torch.sum(_a, axis=-1) * self.mask_for_loss[b_ix]
            _c = self.mask_for_loss[b_ix]
            # tensor
            recovery_rate = torch.sum(_b) / torch.sum(_c)
            recovery_rate = np.float32(recovery_rate.detach().cpu().numpy())
            res.append(recovery_rate)
        # print('recovery_rate:', res)
        return res

    def pairwise_match(self, sample_dict) -> list:
        res = []
        native_seq = [v for k,v in self.protein.items() if k.startswith('seq_chain_')]
        for sample_seq in sample_dict['sample_seq']:
            _scores = 0
            for target, query in zip(native_seq, sample_seq.split('/')):
                alignment = self.align_match(target, query)
                _scores += int(alignment.score)
            res.append(_scores)
        return res

    def align_match(self, seq1, seq2):
        '''
        global alignment
        '''
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 1
        aligner.mismatch_score = 0
        aligner.open_gap_score = 0
        aligner.extend_gap_score = 0
        alignment = aligner.align(seq1, seq2)
        return alignment

    def pairwise_mismatch(self, sample_dict) -> list:
        res = []
        native_seq = [v for k,v in self.protein.items() if k.startswith('seq_chain_')]
        for sample_seq in sample_dict['sample_seq']:
            _scores = 0
            for target, query in zip(native_seq, sample_seq.split('/')):
                alignment = self.align_mismatch(target, query)
                _scores += int(alignment.score)
            res.append(_scores)
        return res

    def align_mismatch(self, seq1, seq2):
        aligner = PairwiseAligner()
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.open_gap_score = 0
        aligner.extend_gap_score = 0
        alignment = aligner.align(seq1, seq2)
        return alignment

    #################################33
    @property
    def X(self):
        self._X[np.isnan(self._X)] = 0.
        X = torch.from_numpy(self._X).to(
            dtype=torch.float32, device=self.device
        )
        return X[:,:,0] if self.conf.ca_only else X

    @property
    def S(self):
        return torch.from_numpy(self._S).to(
            dtype=torch.long,device=self.device
        )

    @property
    def pssm_coef(self):
        return torch.from_numpy(self.pssm_coef_all).to(
            dtype=torch.float32, device=self.device
        )
    
    @property
    def pssm_bias(self):
        return torch.from_numpy(self.pssm_bias_all).to(
            dtype=torch.float32, device=self.device
        )

    @property
    def mask(self):
        _mask = np.isfinite(
            np.sum(self._X,(2,3))).astype(np.float32
        )
        return torch.from_numpy(_mask).to(
            dtype=torch.float32, device=self.device
        )

    @property
    def chain_M(self):
        return torch.from_numpy(self._chain_M).to(
            dtype=torch.float32, device=self.device
        )

    @property
    def chain_M_pos(self):
        return torch.from_numpy(self._chain_M_pos).to(
            dtype=torch.float32, device=self.device
        )

    @property
    def mask_for_loss(self):
        '''
        return torch.Tensor
        '''
        return self.mask * self.chain_M * self.chain_M_pos

    @property
    def tied_pos_list_of_lists_list(self):
        return self._tied_pos_list_of_lists_list

    @property
    def lengths(self):
        #sum of chain seq lengths
        return np.array([len(batch['seq']) for b in self.batch], dtype=np.int32)

    @property
    def omit_AA_mask(self):
        return torch.from_numpy(self._omit_AA_mask).to(
            dtype=torch.float32, device=self.device
        )

    @property
    def tied_beta(self):
        return self._tied_beta
    
    @property
    def bias_by_res_all(self):
        return torch.from_numpy(self._bias_by_res_all).to(
            dtype=torch.float32, device=self.device
        )
    
    @property
    def pssm_log_odds_all(self):
        return torch.from_numpy(self._pssm_log_odds_all).to(
            dtype=torch.float32, device=self.device
        )
    
    @property
    def masked_list_list(self):
        return self._masked_list_list
    
    @property
    def chain_id_list_list(self):
        return self._chain_id_list_list

    @property
    def chain_encoding_all(self):
        return torch.from_numpy(self._chain_encoding_all).to(
            dtype=torch.long, device=self.device
        )
    
    @property
    def visible_list_list(self):
        return self._visible_list_list
    
    @property
    def masked_chain_length_list_list(self):
        return self._masked_chain_length_list_list

    @property
    def residue_idx(self):
        return torch.from_numpy(self._residue_idx).to(
            dtype=torch.long,device=self.device
        )

    @property
    def dihedral_mask(self):
        return torch.from_numpy(self._dihedral_mask).to(
            dtype=torch.float32, device=self.device
        )
    
    @property
    def name(self):
        return self.batch[0]['name']

    @property
    def model_name(self):
        return self.conf.model_name
