"""
    based on https://github.com/ilanschnell/bitarray/blob/master/examples/smallints.py
"""

from pathlib import Path
import re
from bitarray import bitarray

base_seq_types = {
    'dna_bases': {'A': bitarray('00'), 'C': bitarray('01'), 'G': bitarray('10'), 'T': bitarray('11')},
    'rna_bases': {'A': bitarray('00'), 'C': bitarray('01'), 'G': bitarray('10'), 'U': bitarray('11')}
}


class NASeq(object):  # for Nucleic Acid sequences (DNA/RNA) 4 letter alphabets
    """
    sequence type which encodes bases into 2 bit objects
    stored in the largest bit integer arrays for your platform
    """
    __slots__ = ['seq', 'base_seq_type', 'base_seq_type_name', 'metadata']

    def __init__(self, input_seq: str, base_seq_type_name: str, metadata: dict = {}):
        input_seq = input_seq.upper()
        self.base_seq_type_name = base_seq_type_name.lower()
        try:
            self.base_seq_type = base_seq_types[base_seq_type_name]
        except KeyError:
            raise KeyError('The only base_seq_type options are "dna_bases" and "rna_bases"')
        if set(input_seq) - set(self.base_seq_type.keys()):
            raise ValueError('The input_seq contains characters not in the base_seq_type')
        self.seq: bitarray = bitarray()
        self.metadata: dict = metadata
        self.seq.encode(self.base_seq_type, input_seq)

        # number of base pairs
        self.metadata['len'] = int(len(self.seq) / 2)

        # masks the bits at the end to make sure they aren't used for anything.
        # TODO: replace the 32 with the longest int register for this platform

    def __repr__(self):
        letter_seq = ''.join(self.seq.decode(self.base_seq_type))
        return f'NASeq({letter_seq}, {self.base_seq_type_name}'

    def __str__(self):
        return ''.join(self.seq.decode(self.base_seq_type))

    def slice_i(self, i):
        assert 0 <= i < self.metadata['len']
        return slice(self.k * i, self.k * (i + 1))

    def __getitem__(self, i):
        if isinstance(i, slice):
            start = i.start or 0
            stop = i.stop or len(self.seq)
            assert -1 * self.metadata['len'] <= start <= self.metadata['len']
            assert -1 * self.metadata['len'] <= stop <= self.metadata['len']
            if i.step is not None:
                raise NotImplementedError('skips and reverse indexing is not Implemented')
            return ''.join(bitarray(self.seq[start * 2: stop * 2]).decode(self.base_seq_type))

        else:
            if 0 <= i <= self.metadata['len'] - 1:
                return self.seq[i * 2: i * 2 + 2].decode(self.base_seq_type)[0]
            elif -1 * self.metadata['len'] <= i < 0:
                return self.seq[self.metadata['len'] * 2 + (i * 2): self.metadata['len'] * 2 + i * 2 + 2].decode(
                    self.base_seq_type)[0]
            else:
                raise IndexError

    def __setitem__(self, i, v):
        assert 0 <= v < 2 ** self.k
        a = bitarray(endian='little')
        a.fromstring(chr(v))
        self.data[self.slice_i(i)] = a[:self.k]

    def complement(self, i=None):
        if isinstance(i, slice):
            return ''.join(self.seq[i].__xor__(bitarray([1] * len(self.seq[i]))).decode(self.base_seq_type))
        return ''.join(self.seq.__xor__(bitarray([1] * len(self.seq))).decode(self.base_seq_type))

    def __eq__(self, other):
        return self.seq == other.seq

    def __format__(self, format_spec):
        return str(self)

    def __hash__(self):
        return hash(self.seq)

    def __len__(self):
        """Return the length of the sequence, use len(my_seq)."""
        return self.metadata['len']  # Seq API requirement


class FastaData():
    def __init__(self, path):
        self.data = Path(path).read_text()
        self._sequences = self.sequences()
        # self._metadata

    def num_records(self):
        comp = re.compile('^>\w', re.MULTILINE)
        return len(re.findall(comp, self.data))

    def sequences(self):
        seq = ''
        for line in self.data.split('\n'):
            if line.startswith('>'):
                if seq:
                    yield seq
                seq = ''
            else:
                seq += line
        yield seq

    def sequence_lengths(self):
        sequences = []

        return [len(x) for x in self.sequences()]


if __name__ == '__main__':
    pass
    # a = NASeq('acgttgca', 'dna_bases')
    # print(a[2])
    # fd = FastaData('/Users/jlangley/Documents/code/biostuff/datafiles/dna.example.fasta')
    # fd.num_records()
    # fd.sequence_lenths()
    # max(fd.sequence_lengths())
    # min(fd.sequence_lengths())
