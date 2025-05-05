'''
'''
import os
import re
import pickle

class AnalyzeComplex:
    def __init__(self, pickle_file, delta_threshold=0, max_span=30):
        self.load_data(pickle_file)
        self.delta_threshold = delta_threshold
        self.max_span = max_span

    def load_data(self, infile):
        '''
        self.data
        self.outdir
        '''
        self.outdir = os.path.dirname(infile)
        with open(infile, 'rb') as f:
            self.data = pickle.load(f)
        print('number of complex', len(self.data))

    def print_meta(self):
        '''
        print delta SASA
        '''
        for _data in self.data:
            k1, k2 = _data['type'].split('-')
            c1, c2 = _data['chains'][:2]
            print(f'{k1} chain {c1[0]} ~ {k2} chain {c2[0]}')
    
    def print_delta_sasa(self, ix:int):
        '''
        print delta SASA
        '''
        k1, k2 = self.data[ix]['type'].split('-')
        c1, c2 = self.data[ix]['chains'][:2]
        print(f'complex: {k1} chain {c1[0]} ~ {k2} chain {c2[0]}')


        print('---'*10)
        print(f'{k1} chain:')
        delta1 = self.data[ix]['delta_total_sasa']['delta1']
        sig1 = delta1[delta1['value'] > self.delta_threshold]
        print(sig1)

        print('---'*10)
        print(f'{k2} chain:')
        delta2 = self.data[ix]['delta_total_sasa']['delta2']
        sig2 = delta2[delta2['value'] > self.delta_threshold]
        print(sig2)
    
    def retrieve_seq(self):
        pairs = []
        for _data in self.data:
            k1, k2 = _data['type'].split('-')
            # delta df
            delta1 = _data['delta_total_sasa']['delta1']
            delta2 = _data['delta_total_sasa']['delta2']
            # seq
            motif1 = self.parse_seq(delta1)
            motif2 = self.parse_seq(delta2)
            pairs.append({
                k1: motif1,
                k2: motif2,
            })
        # export
        outfile = os.path.join(self.outdir, 'seq_sasa.p')
        with open(outfile, 'wb') as f:
            pickle.dump(pairs, f, protocol=pickle.HIGHEST_PROTOCOL)
        return pairs, outfile


    def parse_seq(self, delta):
        '''
        '''
        sig = delta[delta['value'] > self.delta_threshold]

        # group significant point
        group, curr = [], []
        for ix, res_no in sig['res_no'].items():
            if curr: 
                span = self.cal_span(curr[-1][1], res_no)
                if span <= self.max_span:
                    curr.append((ix, res_no))
                else:
                    group.append(curr)
                    curr = [(ix, res_no),]
            else:
                curr.append((ix, res_no))
        else:
            if curr:
                group.append(curr)

        # convert to seq
        motif = []
        for g in group:
            start, end = g[0][0], g[-1][0]
            sub = delta.loc[start:end,:]
            seq = ''.join(sub['aa'])
            motif.append({
                'seq': seq,
                'sig_res': len(g),
                'start': g[0],
                'end': g[-1],
            })
        return motif

    def cal_span(self, a, b):
        a = int(re.findall(r'\d+', a)[0])
        b = int(re.findall(r'\d+', b)[0])
        return b - a