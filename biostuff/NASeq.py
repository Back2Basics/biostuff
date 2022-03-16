"""
grovers search inverts the function and can find something in a list.
n^(1/2)
"""
import re
from pathlib import Path
from typing import Generator, Tuple

import numpy as np
from more_itertools import chunked

from utilities.compiled_sequences import groups_of_4, check_base_type, seq2base_map

# from toolz import curried as c # check out https://learning.oreilly.com/library/view/elegant-scipy/9781491922927/ch08.html#idm140212611478288

how_to_use_args = "" \
                  "input_seq = dna or rna sequence you want put in using the 4 nucleic acid letters \n" \
                  "base_seq_type = either 'rna_bases' or 'dna_bases' \n"


class NASeq:  # for Nucleic Acid sequences (DNA/RNA) 4 letter alphabets
    """
    sequence type which encodes bases into 2 bit objects
    stored in the largest bit integer arrays for your platform
    """
    # __slots__ = ['seq', 'input_seq', 'seq_len', 'seq_name', 'seq_to_bases', 'bases_to_seq', 'seq_start',
    #              'base_seq_type',
    #              'base_seq_type_name']

    size = 32  # bases per int (64 bit ints on intel based machines)
    # TODO: replace the 32 with the longest int register for this platform
    metadata = dict()
    seq_start = None
    seq_len = None

    def __init__(self, **kwargs):
        """
        this takes a sequence and converts it to a sequence of 2 bit integers
        using the groups_of_4 dictionary to convert the bases to integers
        then stored in a numpy array of ints
        :param kwargs:
        """
        self.seq_to_bases = None
        self.bases_to_seq = None
        self.__dict__.update(kwargs)
        self.input_seq = kwargs.get("input_seq", '').upper()
        self.base_seq_type = kwargs.get("base_seq_type", '').lower()  # "dna_bases" or "rna_bases"
        if not all([self.input_seq, self.base_seq_type]):
            raise ValueError(f"either input_seq or base_seq_type not specified. \n {how_to_use_args}")
        self.seq_len = len(self.input_seq)
        self.check_input_sequence()
        self.set_seq_start()
        self.set_base_seq_type(self.base_seq_type)

        self.input_seq = self.input_seq.encode('ascii')
        self.input_seq_with_leading_zeros = b'A' * self.seq_start + self.input_seq
        # take self.input_seq_with_leading_zeros in groups of 4 letters at a time in chunks

        self.seq = self.seq_array(
            base_seq_type='dna_bases',
            bases_to_seq_or_seq_to_bases="bases_to_seq",
            input_seq=self.input_seq_with_leading_zeros)

        del self.input_seq
        # self.seq = self.convert_bases_to_seq(self.input_seq)

    def seq_array(self,
                  base_seq_type: str = 'dna_bases',
                  bases_to_seq_or_seq_to_bases: str = "bases_to_seq",
                  input_seq: str = ""):
        return np.array(
            [groups_of_4[base_seq_type][bases_to_seq_or_seq_to_bases][input_seq[i:i + 4]] for i in
             range(0, self.seq_len, 4)],
            dtype=np.uint8)

    def set_seq_start(self):
        # checked
        """
        sometimes the seq_start is not 0 and the seq_len is not a multiple of the size
        so this function will set the seq_start to the correct value by
        returning the number of leading zeros needed to make the sequence divisible by size

        :param seq_len: length of the sequence
        :param size: size of the ints in the sequence

        :return: the number of leading 2 bits needed to make the sequence divisible by size modulo size
        """
        if (self.seq_len % self.size) == 0:
            self.seq_start = 0
        else:
            self.seq_start = self.size - (self.seq_len % self.size)

    def check_input_sequence(self):
        # checked
        """
        base_seq_type = "ACGT"
        input_seq = "A"
        """
        if set([ord(x.encode('ascii')) for x in self.input_seq]) - set(
                [ord(x.encode('ascii')) for x in check_base_type[self.base_seq_type]]):
            raise ValueError(
                "The input_seq contains characters not in the right base_type.  e.g. DNA doesn't have 'U' for instance")

    def set_base_seq_type(self, base_seq_type):
        # checked
        try:
            self.bases_to_seq = groups_of_4[self.base_seq_type]['bases_to_seq']
            self.seq_to_bases = groups_of_4[self.base_seq_type]['seq_to_bases']
        except KeyError:
            raise KeyError('The only base_seq_type options are "dna_bases" and "rna_bases"')

    def save_seq(self, filename: Path = Path('sequence.npz')):
        # checked
        np.savez_compressed(filename, self.input_seq_with_leading_zeros)  # saves with a npz extension

    def load_seq(self, filename: Path = Path("sequence.npz")):
        # checked
        with np.load(filename) as data:
            self.input_seq_with_leading_zeros = data['arr_0']
        # self.input_seq_with_leading_zeros = np.load(filename.absolute())

    def __getitem__(self, i):
        if isinstance(i, slice):
            start = i.start or 0
            stop = i.stop or len(self.seq)
            step = i.step or 1
        else:
            if len(self.seq) > 1:
                start = i or 0
                start_int = int(np.floor((start + self.seq_start) / 32))
                start_position = (start + self.seq_start) % 32
                return self.decode_int64(self.seq[start_int:start_int + 1].view())[start_position]
            else:
                start = i or 0
                # start_int = int(np.floor((start + self.seq_start) / 32))
                start_position = (start + self.seq_start) % 32
                return self.decode_int64(self.seq)[start_position]

        start_int = int(np.floor((start + self.seq_start) / 32))
        end_int = int(np.floor((stop + self.seq_start) / 32))
        start_position = (start + self.seq_start) % 32
        end_position = (stop + self.seq_start) % 32
        assert -1 * self.seq_len <= start <= self.seq_len, ValueError(f"value needs to be inside the sequence length")
        assert -1 * self.seq_len <= stop <= self.seq_len, ValueError(f"value needs to be inside the sequence length")
        if step is not 1:
            raise NotImplementedError('skips and reverse indexing is not Implemented')
        return self.convert_seq_to_bases(self.seq[start_int, end_int])[start_position:end_position]

    def print_seq(self, np_array):
        # todo try timing str.maketrans instead of dictionary lookup
        return ''.join([seq2base_map[y] for x in np_array.view(np.uint8) for y in np.base_repr(x, base=4)])

    def __repr__(self):
        # checked
        return f'NASeq({self.base_seq_type} {self.seq_len} bases long)'

    def __str__(self):
        return f'NASeq({self.base_seq_type} {self.seq_len} bases long)'

    def _slice_i(self, i):
        """since each base takes 2 bits everything that needs to be indexed needs to be multiplied by 2
        i = 5 gets converted to slice(10, 12)
        8:10 gets converted to slice(16, 20)
        """
        if isinstance(i, slice):
            return slice(i.start * 2, i.stop * 2)
        else:
            assert 0 <= i < self.seq_len
            return slice(i * 2, i * 2 + 2)

    def __setitem__(self, i, v):
        if i.step != 1:
            raise NotImplemented("step can't be set")
        if not isinstance(i, slice):
            a = self.single_base_to_seq(v)
            self.seq[i.start, i.stop] = a
        else:
            a = self.single_base_to_seq(v)
            self.seq[i] = a

    def seq_complement(self):
        return ~self.seq

    def base_complement(self, i=None):
        print(self.seq_complement())

    def __eq__(self, other):
        return self.seq == other.seq

    def __format__(self, format_spec):
        return str(self)

    def __hash__(self):
        return hash(self.seq)

    def __len__(self):
        """Return the length of the sequence, use len(my_seq)."""
        return self.metadata['len']  # Seq API requirement

    def __contains__(self, seq_fragment):
        lookup_base_seq = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        self.other_fragment_seq = sum([lookup_base_seq[base] * 4 ** count for count, base in enumerate(seq_fragment)])

        # np.left_shift(40, 2)
        # return any(self.seq[2 * i:2 * i + length_a] == a for i in range(len(self.seq)) if i % 2 == 0)

        # def generate_start_codon_locations(self):
        #     """is AUG in self.seq
        #     this can only happen if this is an RNA
        #     """
        #     if self.base_seq_type != base_seq_types["rna_bases"]:
        #         raise ValueError(f"Not an RNA Sequence")
        #     comparison = bitarray().encode(self.base_seq_type, 'AAA')
        #     a = bitarray()
        #     a.encode(self.base_seq_type, 'AUG')
        #     length_a = len(a)
        #     for i in range(self.metadata['len']) - a.metadata['len']:
        #         if (self.seq[
        #             2 * i:2 * i + length_a] ^ a) == comparison:  # TODO: figure out if this is faster than... if x==y for various bit lengths
        #             yield i

    def chunker_for_32_bases(self, seq: str):
        self.seq_len = len(seq)
        prepend_length = 32 - (self.seq_len % 32)
        prepend = b'A' * prepend_length
        for num in it.count(step=32):
            yield it.islice(prepend + seq, num, num + 32)

    def chunker_for_bases(self, seq: str, chunk_size: int) -> Generator[
        Tuple[str, int], str, None]:
        """you can pass in a whole sequence or 32 bases to be broken up into smaller sections"""
        seq_length = len(seq)

        if passed_in_seq_length % chunk_size != 0:
            first_partition = chunk_size - self.seq_start
        else:
            first_partition = chunk_size
        yield seq[0:first_partition] + b''.join(
            [b'A'] * (chunk_size - first_partition))  # f"{seq[0:first_partition]:A>{chunk_size}}"
        if passed_in_seq_length > chunk_size:
            for pos in range(first_partition, passed_in_seq_length, chunk_size):
                yield seq[pos:pos + chunk_size]

    @classmethod
    def little_indian_reverse(cls, seq, chunk_size: int):  # pass in np.uint8 from the np.uint64
        for group in range(0, len(seq), chunk_size):
            seq[group:group + chunk_size] = seq[group:group + chunk_size][::-1]

    def split_into(self, iterable_seq, step=None):
        # checked
        """
        sizes = [5]+ it.repeat(32)
        iterable = b'ACAACACACCAACCCAAACACAC'
        """
        yield from (bytes(x) for x in chunked(iterable_seq, step))

    def convert_bases_to_seq(self, input_seq: bytes):
        """
        turn the
        into a numpy array of int64 types
        """
        print(input_seq)
        val = [y for y in self.split_into(input_seq, 32)]
        # val = np.array([
        #     self.bases_to_seq.__getitem__(y) for x in self.split_into(iterable_seq=input_seq, sizes=32) for y in
        #     self.split_into(iterable_seq=x, sizes=4)])
        return val

    # def convert_bases_to_seq(self, input_seq: str):
    #     val = np.array([
    #         self.bases_to_seq.__getitem__(y) for x in
    #         self.chunker_for_bases(seq=input_seq, chunk_size=32, passed_in_seq_length=self.seq_len) for y in
    #         self.chunker_for_bases(seq=x, chunk_size=4)],
    #         dtype=np.uint8)
    #     self.little_indian_reverse(val, 8)
    #     return val.view(np.uint64)

    @classmethod
    def test_convert_bases_to_seq(cls, input_seq: str, seq_len: int):
        val = np.array([
            cls.bases_to_seq.__getitem__(y) for x in
            cls.chunker_for_bases(seq=input_seq, chunk_size=32, passed_in_seq_length=seq_len) for y in
            cls.chunker_for_bases(seq=x, chunk_size=4)],
            dtype=np.uint8)
        cls.little_indian_reverse(val, 8)
        return val.view(np.uint64)

    def convert_seq_to_bases(self, seq=None):
        # TODO time the bitshifted version to this and see which is faster
        return b''.join((self.seq_to_bases[x] for x in seq.view(dtype=np.uint8)))

    # def has_end_codon(self):
    #     any()


class FastaData():
    def __init__(self, path):
        self.data = Path(path).read_text()
        self._sequences = self.sequences()

    def num_records(self):
        comp = re.compile('^>w')
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
    # todo: incorporate http://hplgit.github.com/bioinf-py/doc/src/data/genetic_code.tsv

    # a = NASeq('acgttgca', 'dna_bases')
    # print(a[2])
    # fd = FastaData('/Users/jlangley/Documents/code/biostuff/datafiles/dna.example.fasta')
    # fd.num_records()
    # fd.sequence_lenths()
    # max(fd.sequence_lengths())
    # min(fd.sequence_lengths())
