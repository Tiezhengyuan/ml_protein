import os
import numpy as np
import pandas as pd


class RetrievePdb:

    def __init__(self, conf):
        self.conf = conf

    def homomer_tied_positions(self) -> dict:
        '''
        title Helper functions
        '''
        my_dict = {}
        for result in self.conf['pdb_dict_list']:
             #A, B, C, ...
            all_chain_list = sorted([
                item[-1:] for item in list(result) if item[:9] == 'seq_chain'
            ])
            tied_positions_list = []
            chain_length = len(result[f"seq_chain_{all_chain_list[0]}"])
            for i in range(1,chain_length+1):
                temp_dict = {}
                for j, chain in enumerate(all_chain_list):
                    temp_dict[chain] = [i] #needs to be a list
                tied_positions_list.append(temp_dict)
            my_dict[result['name']] = tied_positions_list
        return my_dict