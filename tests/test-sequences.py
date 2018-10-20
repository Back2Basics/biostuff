from Python_for_Bioinformatics_final import NASeq
import unittest as ut


class ReadOnlyTests(ut.TestCase):
    def setUp(self):
        self.dna = NASeq('ACGTTGCA', 'dna_bases')

    def test__getitem__(self):
        assert self.dna[0] == 'A'
        assert self.dna[1] == 'C'
        assert self.dna[2] == 'G'
        assert self.dna[-1] == 'A'
        assert self.dna[-3] == 'G'
        # assert self.dna[-1:-3:-1] == 'CG'

    def test_slice(self):
        assert self.dna[2:5] == 'GTT'
        assert self.dna[:2] == 'AC'
        assert self.dna[1:3] == 'CG'

    def test_reverse_complement(self):
        self.assertEqual(self.dna.complement(), 'TGCAACGT')


# class WriteableTests(ut.TestCase):
#     def setUp(self):
#         self.dna = NASeq('ACGTTGCA', 'dna_bases')
#
#     def test__setitem__individual(self):
#         self.dna[0] = 'T'
#         assert self.dna[0]=='T'
#
#     def test__setitem__slice(self):
#         self.dna[1:3] = 'AT'
#         assert self.dna[1:3] == 'AT'
#
