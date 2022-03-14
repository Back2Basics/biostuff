import itertools as it

import numpy as np

dna_bases = {'A', 'C', 'G', 'T'}
rna_bases = {'A', 'C', 'G', 'U'}
check_base_type = {"dna_bases": dna_bases, "rna_bases": rna_bases}
seq2base_map = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
# dna_str_trans = str.maketrans('0123', b'ACGT')
# rna_str_trans = str.maketrans('0123', b'ACGU')


def build_numpy_array_with_ones(shape):
    """return a numpy array with 1's in the shape of shape"""


# create 8 bit lookup tables
def build_seq_to_bases(base_seq_type):
    return np.array([''.join(x) for x in it.product(base_seq_type, repeat=4)], dtype="|S4")


def build_bases_to_seq_dict(base_seq_type):
    return {x: count for count, x in enumerate(build_seq_to_bases(base_seq_type))}


four_dna_bases_to_seq_combos = build_bases_to_seq_dict(dna_bases)
four_rna_bases_to_seq_combos = build_bases_to_seq_dict(rna_bases)

np_seq_to_four_dna_bases_combos = build_seq_to_bases(dna_bases)
np_seq_to_four_rna_bases_combos = build_seq_to_bases(rna_bases)

groups_of_4 = {
    'dna_bases':
        {'bases_to_seq': four_dna_bases_to_seq_combos,
         'seq_to_bases': np_seq_to_four_dna_bases_combos},
    'rna_bases':
        {'bases_to_seq': four_rna_bases_to_seq_combos,
         'seq_to_bases': np_seq_to_four_rna_bases_combos}
}
