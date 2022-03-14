import pytest

from biostuff.NASeq import *


def test_NASeq__repr__():
    seq = NASeq(input_seq='A', base_seq_type='dna_bases')
    assert f'NASeq({seq.base_seq_type} {seq.seq_len} bases long)' == seq.__repr__()


def test_set_base_seq_type():
    seq = NASeq(input_seq='A', base_seq_type='dna_bases')
    seq.set_base_seq_type('dna_bases')
    assert b'AAAA' in seq.bases_to_seq
    # assert 170 in seq.seq_to_bases


def test_set_seq_start_test():
    seq = NASeq(input_seq='A', base_seq_type='dna_bases')
    assert seq.seq_start == 31
    seq = NASeq(input_seq="A"*20, base_seq_type='dna_bases')
    assert seq.seq_start == 12
    seq = NASeq(input_seq='A'*33, base_seq_type='dna_bases')
    assert seq.seq_start == 31
    seq = NASeq(input_seq='A'*33, base_seq_type='dna_bases')
    assert not seq.seq_start == 30
    seq = NASeq(input_seq='A'*65, base_seq_type='dna_bases', size=64)
    assert seq.seq_start == 63

def test_save_seq_and_load_seq():
    seq = NASeq(input_seq='AGA', base_seq_type='dna_bases')
    filename = Path("test_save_seq_and_load_seq.npz")
    if filename.exists():
        filename.unlink()

    seq.save_seq(filename=filename)
    assert filename.exists()
    seq.load_seq(filename)
    assert seq.input_seq == b'AGA'
    filename.unlink()

def test_required_NASeq_elements():
    with pytest.raises(ValueError):
        seq = NASeq(input_seq='A')
    with pytest.raises(ValueError):
        seq = NASeq(dna_seq_bases='dna_bases')


def test_NASeq_stuff():
    seq = NASeq(input_seq='CA', base_seq_type="dna_bases")
    assert seq.input_seq == b'CA'


def test_NASeq_leading_zeroes():
    seq = NASeq(input_seq='CA', base_seq_type="dna_bases")
    assert seq.input_seq_with_leading_zeros == b'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAACA'
    seq = NASeq(input_seq='A', base_seq_type="dna_bases")
    assert seq.input_seq_with_leading_zeros == b'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
    seq = NASeq(input_seq='C' * 32, base_seq_type="dna_bases")
    assert seq.input_seq_with_leading_zeros == b'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
    seq = NASeq(input_seq='C' * 33, base_seq_type="dna_bases")
    assert seq.input_seq_with_leading_zeros == b'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
    seq = NASeq(input_seq='ACGTTGCT', base_seq_type="dna_bases")
    assert seq.input_seq_with_leading_zeros == b'AAAAAAAAAAAAAAAAAAAAAAAAACGTTGCT'
    seq = NASeq(input_seq='ACGTTGCT', base_seq_type="dna_bases")
    assert seq.input_seq_with_leading_zeros != b'ACGTTGCT'


def test_NASeq_check_input_sequences():
    with pytest.raises(ValueError):
        c = NASeq(input_seq='UA', base_seq_type='dna_bases')

    with pytest.raises(ValueError):
        c = NASeq(input_seq='A-', base_seq_type='dna_bases')

    b = NASeq(input_seq='CA', base_seq_type='dna_bases')
    # if set([ord(x.encode('ascii')) for x in self.input_seq]) - set([ord(x.encode('ascii')) for x in base_type[self.base_seq_type]]):
    #     raise ValueError("The input_seq contains characters not in the right base_type.  e.g. DNA doesn't have 'U' for instance")


# def test_NASeq_convert_bases_to_seq():
#     c = NASeq(input_seq='CA', base_seq_type='dna_bases')
#     assert 'C' == c.convert_bases_to_seq(b'AAAAAAAAAAAAAAAAAAAAAAAAACGTTGCT')
#     c = NASeq(input_seq='C'*33, base_seq_type='dna_bases')
#     assert 'C' == c.convert_bases_to_seq(b'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC')
#

def test_NASeq_split_into():
    seq = NASeq(input_seq='CAAAAAAAAAAAAAAAAAAAAAATTAAAAAAACA', base_seq_type="dna_bases")
    assert list(seq.split_into(b'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAATTAAAAAAACA', 32)) == [
        b'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAACA', b'AAAAAAAAAAAAAAAAAAAAATTAAAAAAACA']

    # @mark.parametrize('test_input, expected', [('A', b'A')])
    # def test_input_type_test(test_input, expected):
    #     b = NASeq.check_input_sequence_test(base_seq_type='ACGT', input_seq=test_input)
    #     assert b == expected

    # @mark.parametrize(
    #     "test_input,expected",
    #     [('A', b'A'), ('AC', b'AC')])
    # def test_input_type(test_input, expected):
    #     a = NASeq(input_seq=test_input, base_seq_type='dna_bases')
    #     b = a.check_input_sequence()
    #     assert b == expected
    #
    #
    # @mark.parametrize("input_test,expected", [('A', b'A')])
    # def test_convert_bases_to_seq(input_test, expected):
    #     a = NASeq.test_convert_bases_to_seq(input_test, len(input_test))
    #     assert a == expected
    #
    #
    # @mark.parametrize('input_test,expected', [(b'A', b'AAAA'),
    #                                                  (b'CA', b'CAAA'),
    #                                                  (b'AC', b'ACAA'),
    #                                                  (b'CAC', b'CACA'),
    #                                                  (b'GT', b'GTAA'),
    #                                                  (b'GTTT', b'GTTT'),
    #                                                  (b'GTTTA', b'GTTTAAAA'),
    #                                                  (b'GTTTCACCAG', b'GTTTCACCAGAA'),
    #                                                  ])
    # def test_chunker_for_bases(input_test, expected):
    #     a = list(NASeq.test_chunker_for_bases(seq=input_test, chunk_size=4, seq_start=4 - len(input_test)))
    #
    #     assert a == list(bytes(x) for x in more_itertools.grouper(expected, 4))
    #
    #
    # @mark.parametrize(
    #     "test_input,expected",
    #     [(('A', 32), [''.join(['AAAA'] * 8)])])
    # def test_chunker_for_bases_pass(test_input, expected):
    #     a = NASeq(test_input[0], base_seq_type="dna_bases")
    #     assert list(a.chunker_for_bases(
    #         seq=test_input[0],
    #         chunk_size=test_input[1],
    #         passed_in_seq_length=len(test_input[0]))) == expected
    #
    #
    # @mark.parametrize("seq_len,size,expected", [(1, 32, 31), (31, 32, 1), (32, 32, 0)])
    # def test_set_seq_start_test(seq_len, size, expected):
    #     output = NASeq.set_seq_start_test(seq_len, size)
    #     assert expected == output
    #
    #
    # @mark.parametrize("test_input,expected",
    #                          [(string.ascii_letters,
    #                            ValueError),
    #                           (None, TypeError)
    #                           ])
    # def test_init_failures1(test_input, expected):
    #     with pytest.raises(expected):
    #         a = NASeq(test_input, base_seq_type="dna_bases")
    #
    #
    # @mark.parametrize("test_input,expected",
    #                          [("something_wrong",
    #                            KeyError)
    #                           ])
    # def test_init_failures2(test_input, expected):
    #     with pytest.raises(expected):
    #         a = NASeq("AAA", base_seq_type="test_input")
    #
    #
    # @mark.parametrize(
    #     "test_input,expected",
    #     [
    #         (("AAAAAAAGCGTTTGCCAACCAAAGCGTTTGCCAACCAAA", 32), 25),
    #         (("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGCGTTTGCCAACCAAAGCGTTTGCCAACCAAA",
    #           32), 0),
    #         (("AAAAAAAAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT", 32), 25),
    #         (("AAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", 32), 25),
    #         (("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTGCA", 32), 0),
    #     ])
    # def test_start_position(test_input, expected):
    #     a = NASeq(test_input[0], base_seq_type="dna_bases")
    #     assert a.seq_start == expected
    #
    #
    # @mark.parametrize(
    #     "test_input,expected",
    #     [(("A", 32), [0]),
    #      (("AGCGTTTGCCAACCAAAGCGTTTGCCAACCAA", 32), [2809771522707574864]),
    #
    #      (("AAAAAAAGCGTTTGCCAACCAAAGCGTTTGCCAACCAAA", 32), [0, 11239086090830299456]),
    #      (("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGCGTTTGCCAACCAAAGCGTTTGCCAACCAAA",
    #        32), [18446744073709551615, 18446744073709551615, 11239086090830299456]),
    #      (("AAAAAAAAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT", 32), [0, 24113385219402495]),
    #      (("AAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", 32), [0, 18446744073709551615]),
    #      (("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTGCA", 32), [18446744073709551588]),
    #      ])
    # def test_convert_32_bases(test_input, expected):
    #     """ tests to see if the little indian order is removed"""
    #     a = NASeq(test_input[0], base_seq_type="dna_bases")
    #     np.testing.assert_equal(a.seq, np.array(expected, dtype=np.uint64))
    #
    #
    # @mark.parametrize(
    #     "test_input,expected",
    #     [('G', np.array([144115188075855872], dtype=np.uint64)),
    #      ])
    # def test_little_indian_reverse(test_input, expected):
    #     a = NASeq(test_input, base_seq_type="dna_bases")
    #     NASeq.little_indian_reverse(a.seq.view(np.uint8), chunk_size=8)
    #     np.testing.assert_equal(expected, a.seq)
    #
    #
    # @mark.parametrize(
    #     "test_input, expected", [
    #         (0, b'A'),
    #         (1, b'C'),
    #         (2, b'G'),
    #         (3, b'T')]
    # )
    # def test___getitem__(test_input, expected):
    #     a = NASeq('ACGT', base_seq_type="dna_bases")
    #     assert a[test_input] == expected
    #
    #
    # @mark.parametrize(
    #     "test_input, expected", [
    #         (32, b'A'),
    #         (33, b'C'),
    #         (34, b'G'),
    #         (35, b'T')]
    # )
    # def test___getitem__2(test_input, expected):
    #     a = NASeq('ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT', base_seq_type="dna_bases")
    #     assert a[test_input] == expected
    #
    #
    # @mark.parametrize(
    #     "test_input, expected", [
    #         (0, slice(0, 2, None)),
    #         (1, slice(2, 4, None)),
    #         (slice(1, 3), slice(2, 6, None)),
    #         (slice(0, 2), slice(0, 4, None))
    #     ])
    # def test_slice_i(test_input, expected):
    #     a = NASeq('ACGT', base_seq_type="dna_bases")
    #     assert a._slice_i(test_input) == expected
    #
    # #         assert self.dna[1] == 'C'
    # #         assert self.dna[2] == 'G'
    # #         assert self.dna[-1] == 'A'
    # #         assert self.dna[-3] == 'G'
    # #         # assert self.dna[-1:-3:-1] == 'CG'
    # #
    # #     def test_slice(self):
    # #         assert self.dna[2:5] == 'GTT'
    # #         assert self.dna[:2] == 'AC'
    # #         assert self.dna[1:3] == 'CG'
    # #
    # #     def test_reverse_complement(self):
    # #         self.assertEqual(self.dna.complement(), 'TGCAACGT')
    # #
    # #
    # # class WriteableTests():
    # #     def setUp(self):
    # #         self.dna = NASeq('ACGTTGCA', 'dna_bases')
    # #
    # #     def test__setitem__individual(self):
    # #         self.dna[0] = 'T'
    # #         assert self.dna[0] == 'T'
    # #
    # #     def test__setitem__slice(self):
    # #         self.dna[1:3] = 'AT'
    # #         assert self.dna[1:3] == 'AT'
    # #
    # #
    # # class ContainsSeq():
    # #     def setUp(self):
    # #         self.dna = NASeq('ATCTTGCA', 'dna_bases')
    # #
    # #     def test__contains__start_codon(self):
    # #         assert self.dna.__contains__('ATC')
    # #         assert self.dna.__contains__('GCA')
    # #         assert self.dna.__contains__('TGC')
    # #         assert self.dna.__contains__('A')
    # #         assert self.dna.__contains__('C')
    # #         assert self.dna.__contains__('T')
    # #         assert self.dna.__contains__('G')
    # #         assert self.dna.__contains__('ATCTTGCA')
    # #         assert ~self.dna.__contains__('TCA')
    # #
    # #
    # # class StartCodonSeq():
    # #     def setUp(self):
    # #         self.rna = NASeq('AUCUUGCA', 'rna_bases')
    # #
    # #     def test__contains__start_codon(self):
    # #         self.assertEqual(list(self.rna.generate_start_codon_locations()), 10000)
    # #
    # #
    # # class StartCodonSeqErrors():
    # #     def setUp(self):
    # #         self.dna = NASeq('ATCTTGCA', 'dna_bases')
    # #
    # #     def test__contains__start_codon(self):
    # #         with self.assertRaises(ValueError):
    # #             len(list(self.dna.generate_start_codon_locations())) >= 1
