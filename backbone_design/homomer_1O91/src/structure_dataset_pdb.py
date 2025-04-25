'''
build dataset
'''
class StructureDatasetPDB():
    alphabet = 'ACDEFGHIKLMNPQRSTVWYX-'

    def __init__(self, pdb_dict_list, truncate=None, max_length=100):
        alphabet_set = set([a for a in self.alphabet])
        discard_count = {
            'bad_chars': 0,
            'too_long': 0,
            'bad_seq_length': 0
        }

        self.data = []
        for i, entry in enumerate(pdb_dict_list):
            seq = entry['seq']
            name = entry['name']

            bad_chars = set([s for s in seq]).difference(alphabet_set)
            if len(bad_chars) == 0:
                if len(entry['seq']) <= max_length:
                    self.data.append(entry)
                else:
                    discard_count['too_long'] += 1
            else:
                discard_count['bad_chars'] += 1

            # Truncate early
            if truncate is not None and len(self.data) == truncate:
                return

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx]
