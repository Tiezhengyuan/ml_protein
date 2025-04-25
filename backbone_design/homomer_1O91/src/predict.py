import numpy as np
import torch
import torch.nn as nn

from retrieve_pdb import RetrievePdb

class Predict:

    def __init__(self, model, featurizer, conf=None):
        self.model = model
        self.ft = featurizer
        self.conf = conf
        
        #
        self.randn = torch.randn(self.ft.chain_M.shape, device=self.ft.X.device)
        if self.conf:
            self.omit_AAs_np = np.array(
                [AA in self.conf.omit_AAs for AA in self.conf.alphabet]
            ).astype(np.float32)
            self.bias_AAs_np = np.zeros(len(self.conf.alphabet))
            #1.0 for true, 0.0 for false
            self.pssm_log_odds_mask = (
                self.ft.pssm_log_odds_all > self.conf.pssm_threshold
            ).float()
        
        # objects detemined by model execution
        self.sample_dict = None
        self.log_probs = None

    def __call__(self, temperature:int):
        # model sampling: S, probs, decoding_order
        sample_dict = self.cal_sample(temperature)
        sample_dict['temperature'] = [temperature,] * self.conf.batch_size

        # predict probs: shape=(batch_size,393,21)
        sample_dict['log_probs'] = self.cal_sample_probs()
        # scores: (batch_size,)
        sample_dict['score'] = self.cal_scores()
        return sample_dict


    def cal_sample(self, temperature):
        rp = RetrievePdb(self.conf)
        tied_positions_dict = rp.homomer_tied_positions() if self.conf.homomer else None

        if tied_positions_dict == None:
            return self.basic_sample(temperature)
        return self.tied_sample(temperature)


    def basic_sample(self, temperature):
        self.sample_dict = self.model.sample(
            X = self.ft.X,
            randn = self.randn,
            S_true = self.ft.S,
            chain_mask = self.ft.chain_M,
            chain_encoding_all = self.ft.chain_encoding_all,
            residue_idx = self.ft.residue_idx,
            mask = self.ft.mask,
            temperature = temperature, 
            omit_AAs_np = self.omit_AAs_np,
            bias_AAs_np = self.bias_AAs_np, 
            chain_M_pos = self.ft.chain_M_pos,
            omit_AA_mask = self.ft.omit_AA_mask, 
            pssm_coef = self.ft.pssm_coef,
            pssm_bias = self.ft.pssm_bias,
            pssm_multi = self.conf.pssm_multi,
            pssm_log_odds_flag = bool(self.conf.pssm_log_odds_flag), 
            pssm_log_odds_mask = self.pssm_log_odds_mask, 
            pssm_bias_flag = bool(self.conf.pssm_bias_flag), 
            bias_by_res = self.ft.bias_by_res_all
        )
        return self.sample_dict
    
    def tied_sample(self, temperature):
        self.sample_dict = self.model.tied_sample(
            X = self.ft.X,
            randn = self.randn,
            S_true = self.ft.S,
            chain_mask = self.ft.chain_M,
            chain_encoding_all = self.ft.chain_encoding_all, 
            residue_idx = self.ft.residue_idx,
            mask = self.ft.mask,
            temperature = temperature,
            omit_AAs_np = self.omit_AAs_np, 
            bias_AAs_np = self.bias_AAs_np,
            chain_M_pos = self.ft.chain_M_pos, 
            omit_AA_mask = self.ft.omit_AA_mask,
            pssm_coef = self.ft.pssm_coef, 
            pssm_bias = self.ft.pssm_bias,
            pssm_multi = self.conf.pssm_multi, 
            pssm_log_odds_flag = bool(self.conf.pssm_log_odds_flag),
            pssm_log_odds_mask = self.pssm_log_odds_mask, 
            pssm_bias_flag = bool(self.conf.pssm_bias_flag),
            tied_pos = self.ft.tied_pos_list_of_lists_list[0], 
            tied_beta = self.ft.tied_beta, 
            bias_by_res = self.ft.bias_by_res_all
        )
        return self.sample_dict

    def cal_probs(self):
        self.log_probs = self.model(
            X = self.ft.X,
            S = self.ft.S,
            mask = self.ft.mask,
            chain_M = self.ft.chain_M * self.ft.chain_M_pos,
            residue_idx =self.ft.residue_idx,
            chain_encoding_all = self.ft.chain_encoding_all,
            randn = self.randn
        )
        return self.log_probs

    def cal_sample_probs(self):
        self.log_probs = self.model(
            X = self.ft.X,
            S = self.sample_dict["S"],
            mask = self.ft.mask,
            chain_M = self.ft.chain_M * self.ft.chain_M_pos,
            residue_idx = self.ft.residue_idx, 
            chain_encoding_all = self.ft.chain_encoding_all,
            randn = self.randn, 
            use_input_decoding_order = True, 
            decoding_order = self.sample_dict["decoding_order"]
        )
        return self.log_probs

    def cal_scores(self):
        """
        Negative log probabilities
        """
        S = self.sample_dict["S"]
        mask = self.ft.mask_for_loss
        # 
        criterion = torch.nn.NLLLoss(reduction='none')
        loss = criterion(
            self.log_probs.contiguous().view(-1, self.log_probs.size(-1)),
            S.contiguous().view(-1)
        ).view(S.size())
        scores = torch.sum(loss * mask, dim=-1) / torch.sum(mask, dim=-1)
        return scores.cpu().data.numpy()

