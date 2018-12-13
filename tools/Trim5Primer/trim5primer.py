#!/usr/bin/env python
import sys, gzip, re, argparse

from gzip                import GzipFile
from Bio                 import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from itertools           import izip, chain
from textwrap            import dedent
from math                import ceil

__version__ = '1.0.0'

class BetterFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

def parseCommandArguments(*args):
    parser = argparse.ArgumentParser(
        description = dedent('''\
            Clips 5' primers from FastQ paired-end reads

            Primers may be named with a trailing /1 or /2 to apply only to
            forward or reverse reads, respectively.  Reverse complements must
            be explicitly provided if desired.  Primers are compared to the
            first X bases of each read, where X is the length of each primer.
            Walking the primer towards the 3' end of the read is not performed.
            The longest primer to match a read is the one clipped.  Unknown
            bases in reads (designated by N) are always treated as a mismatch.

            Local pairwise alignment is not performed, so indels within the
            primer region of the read are not detected.  If the indel appears
            early in the primer region, it'll likely produce unclipped reads
            due to a high mismatch count.

            Specify file suffixes of .gz to read/write gzipped FastQ.
        '''),
        epilog = dedent('''\
            Trim5Primer version %s
            Copyright 2015 by Thomas Sibley <trsibley@uw.edu>
            Mullins Lab - Microbiology
            University of Washington
        ''' % __version__),
        formatter_class = BetterFormatter)

    for opt in ('fwd.fq', 'forward'), ('rev.fq', 'reverse'):
        parser.add_argument("fastqIn",
            action   = "append",
            help     = "Input FastQ file of %s reads" % opt[1],
            metavar  = opt[0])
    for opt in ('fwd-clipped.fq', 'forward'), ('rev-clipped.fq', 'reverse'):
        parser.add_argument("fastqOut",
            action   = "append",
            help     = "Output FastQ file of clipped %s reads" % opt[1],
            metavar  = opt[0])

    parser.add_argument("-p", "--primers",
        dest     = "fasta",
        required = True,
        help     = "FASTA file of 5' primers",
        metavar  = "clip.fa")
    parser.add_argument("-m", "--mismatches",
        dest     = "mismatches",
        type     = int,
        default  = 3,
        help     = "Number of mismatches allowed when comparing a primer to a read",
        metavar  = "M")
    parser.add_argument("--report-every",
        dest     = 'reportEvery',
        type     = int,
        default  = 100000,
        help     = "Report progress every N read pairs",
        metavar  = 'N')
    return vars(parser.parse_args(*args))

def packSeq(seq, until):
    if until == -1:
        until = len(seq)
    base = {
        'N': 0x0,
        'A': 0x1,
        'T': 0x2,
        'C': 0x4,
        'G': 0x8,
    }
    packed = 0
    for char in seq[0:until]:
        packed = (packed << 4) | base[char]
    return packed

def mismatchCount(bit_int):
    count = 0
    while (bit_int):
        count += (bit_int & 1)
        bit_int >>= 1
    # Each single-base mismatch is two bits XOR'd, except for Ns which are all
    # zeros and will only cause 1 bit flipped.  This produces 0.5 when
    # dividing, which we round up.
    return int(ceil(count / 2.0))

def openMaybeGzip(file, mode):
    f = open(file, mode)
    if re.search('\.gz$', file):
        f = GzipFile(fileobj = f)
    return f

def readPrimers(fasta):
    primers = []
    for primer in SeqIO.parse(open(fasta, "rU"), "fasta"):
        direction = re.search('/([12])$', primer.id)
        direction = [int(direction.group(1))] if direction else [1, 2]
        primers.append(
            (str(primer.seq), packSeq(str(primer.seq), -1), [d - 1 for d in direction]))
    return primers

class Trim5Primer:
    def __init__(self, **args):
        self.primers    = readPrimers(args['fasta'])
        self.mismatches = args['mismatches']

        self.input  = [ openMaybeGzip(i, "r") for i in args['fastqIn']  ]
        self.output = [ openMaybeGzip(o, "w") for o in args['fastqOut'] ]

        (self.r1, self.r2) = [ FastqGeneralIterator(i) for i in self.input ]

        self.count = {
            'pairs':   0,
            'clipped': [0, 0],
        }
        self.reportEvery = args['reportEvery']

    def run(self):
        for reads in izip(self.r1, self.r2):
            self.count['pairs'] += 1

            # Compare each primer to the start (5' end) of each read, recording the
            # maximum offset from 0 that a primer sequence matched within the provided
            # mismatch threshold.  Output each read starting from after the primer, if
            # any.
            for (id, seq, qual), readNum in izip(reads, [0, 1]):
                clipOffset = 0
                for (primer, packedPrimer, direction) in self.primers:
                    if (readNum not in direction) or (len(seq) < len(primer)):
                        continue
                    packedSeq = packSeq(seq, len(primer))
                    if mismatchCount(packedPrimer ^ packedSeq) <= self.mismatches:
                        clipOffset = max(clipOffset, len(primer))
                self.output[readNum].write("@%s\n%s\n+\n%s\n" % (id, seq[clipOffset:], qual[clipOffset:]))
                if clipOffset > 0:
                    self.count['clipped'][readNum] += 1

            if self.count['pairs'] % self.reportEvery == 0:
                print >> sys.stderr, "Processed %d read pairs" % self.count['pairs']

        for handle in chain(self.input, self.output):
            handle.close()

    def reportStats(self):
        print "%i reads (%i pairs) processed" % (self.count['pairs'] * 2, self.count['pairs'])
        if self.count['pairs'] > 0:
            print "%i (%0.1f%%) reads clipped" % (sum(self.count['clipped']), sum(self.count['clipped']) / float(self.count['pairs'] * 2) * 100)
            for (i, direction) in (0, 'forward'), (1, 'reverse'):
                print "%i (%0.1f%%) %s reads clipped" % (self.count['clipped'][i], self.count['clipped'][i] / float(self.count['pairs']) * 100, direction)

if __name__ == '__main__':
    clipper = Trim5Primer( **parseCommandArguments() )
    clipper.run()
    clipper.reportStats()
