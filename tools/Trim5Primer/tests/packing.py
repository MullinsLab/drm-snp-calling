import unittest
from trim5primer import packSeq, mismatchCount

class TestSequencePacking(unittest.TestCase):
    def test_packSeq(self):
        seq = "ATCG"

        # -1 works as a "whole length" flag
        self.assertEqual(
            packSeq(seq, -1),
            packSeq(seq, len(seq)))

        # Max length works
        self.assertEqual(
            packSeq("AT", -1),
            packSeq("ATCG", 2))

        # Packing is working as expected for bitfields
        self.assertEqual(packSeq(seq, -1), 0b0001001001001000)

        # Ns at various positions
        self.assertEqual(packSeq("NAA", -1), 0b00010001)
        self.assertEqual(packSeq("ANA", -1), 0b000100000001)
        self.assertEqual(packSeq("AAN", -1), 0b000100010000)

    def test_mismatchCount(self):
        self.assertMismatchCount("ATCG", "ATCG", 0)
        self.assertMismatchCount("ATCG", "ATTG", 1)
        self.assertMismatchCount("ATCG", "ATTT", 2)
        self.assertMismatchCount("ATCG", "TAGC", 4)
        self.assertMismatchCount("ATCG", "NTCG", 1)
        self.assertMismatchCount("ATCG", "ANCG", 1)
        self.assertMismatchCount("ATCG", "ATNG", 1)
        self.assertMismatchCount("ATCG", "ATCN", 1)

    def assertMismatchCount(self, seqA, seqB, count):
        self.assertEqual(
            mismatchCount(packSeq(seqA, -1)
                        ^ packSeq(seqB, -1)), count)

if __name__ == '__main__':
    unittest.main()
