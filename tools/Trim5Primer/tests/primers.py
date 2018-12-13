import unittest
from trim5primer import readPrimers, packSeq
from tempfile    import NamedTemporaryFile
from textwrap    import dedent

class TestPrimerParse(unittest.TestCase):
    def setUp(self):
        self.fasta = NamedTemporaryFile()
        print >> self.fasta, dedent('''\
            >one/1
            ATCG
            >two/2
            TAGC
            >three
            AT
            TA
        ''')
        self.fasta.flush()

    def test_readPrimers(self):
        primers = readPrimers(self.fasta.name)
        self.assertEqual(primers, [
            ('ATCG', packSeq('ATCG', -1), [0]),
            ('TAGC', packSeq('TAGC', -1), [1]),
            ('ATTA', packSeq('ATTA', -1), [0,1]),
        ])

if __name__ == '__main__':
    unittest.main()
