# Inputs

## Sequencing run for a sample/patient

For each sequencing run, the pipeline expects two raw FastQ files of Illumina
paired-end reads, named like so:

    <...>_R1.fastq.gz       Raw forward reads
    <...>_R2.fastq.gz       Raw reverse reads

These must be available in the reads/ directory, but may be symlinks to
anywhere else.  This lets you preserve the original file names as you wish and
store them wherever you want by simply creating symlinks to use them in the
pipeline.  You can symlink files like so:

    ln -s <existing target> <new symlink>
    ln -s raw/example_A10_0001_R1_ABC.fastq.gz example_R1.fastq.gz

`ln` stands for "link" and `-s` means "create a symbolic (in name only) link".

The slugs <...> in the example fastq filenames above may be whatever you want,
and they will be used as a common prefix for the files which the pipeline will
create.  Runs should have unique prefixes, commonly the patient/sample ID plus
some run ID.

## References and codons of interest

Each sample/patient run must have a reference listed in the
`amplicon-ref-by-run.csv` file in the top-level directory.  The format is:

    run,reference,tagged,comment

The run column must contain the run's name, i.e. the name before `_R1.fastq.gz`
and `_R2.fastq.gz`.

The tagged column should be "yes" or "no" to indicate if the run's data uses
template tags.  The comment may be a short bit of text to note things.

The reference column should contain a reference name.  The reference name must
correspond to the following files:

    refs/<refname>.fasta
    refs/<refname>.gb
    refs/<refname>.bed

The fasta file should contain a single sequence with the same name as
<refname>.  The GenBank (.gb) file should contain a full GenBank entry with
annotations and sequence data limited to the same region as the fasta.  For
example, see refs/AF033819_pol.gb.  The GenBank file is used to build a minimal
SnpEff database.  Kenya references are currently special in that the pipeline
knows how to build a GenBank file for them (assuming nothing changes).

The BED file (.bed) should be a standard UCSC BED format containing lines
defining the codons of interest.  Note that intervals are 0-based, half-open,
meaning that the first three bases in the reference AF033819_pol are
represented as:

    AF033819_pol    0   3   some_codon

The fourth BED column is used as the DRM/codon name in output files.

# Usage

The pipeline uses a Makefile and the command `make` to handle the various
steps.  You run the pipeline, or parts of it, by asking for a specific output
file.  `make` will figure out what inputs it needs and run any steps required.
You can ask for multiple outputs and `make` will produce them in turn.

The final outputs are results/<...>.tsv and summary/<...>.tsv, summary tables of
the called SNPs within a DRM codon and some annotations about each SNP.  To
produce it for the input files 37858_R1.fastq.gz and 37858_R2.fastq.gz,
you run:

    make summary/37858.tsv

This will create results/37858.tsv because the summary tsv depends on it.

The first time you run a new sample, make will need to first produce all the
intermediate files necessary to finally produce the tabular SNP calls for your
FastQ inputs.  Subsequent runs will only produce a new output if one of the
input files has changed.  You can force make to forget about past runs and
output by running:

    make clean

but be aware that this will delete your results!  Luckily, you can regenerate
them easily.

All intermediate files between the raw FastQ files and the final .tsv are
stored in subdirectories of the cache/ directory.  You're free to inspect or
modify these manually.

    cache/
        ref/                Indexed copies of the reference fasta and BED files for each sample
        qc/                 FastQ reads filtered by sickle
        mapped/             mapped sam/bam by bwa mem
        drm-codons/         bam with only mapped reads intersecting a drm codon
        lofreq/             called SNPs via LoFreq*
        haplotypes/         observed haplotypes based on SNPs via Freebayes
        snpeff/             annotated haplotypes/SNPs via SnpEff
        snp-codons/         further annotated haplotypes/SNPs with DRM names
    results/                TSV generated from the snp-codons/ VCFs
    summary/                pared-down version of results TSV

If you want one of the intermediate files without generating the final .tsv,
you can just ask `make` for it.  For example, if you wanted the initial mapping
by bwa mem:

    make cache/mapped/37858.bam

Steps will be re-run whenever the input files are newer than the output files
(or the output files don't exist yet).

## FastQC

You can ask for a FastQC report, as HTML, for a single fastq input like so:

    make 37858_R1_fastqc.html

or a set of fastq files:

    make 37858_R1_fastqc.html 37858_2_fastqc.html

The HTML file is all you need to view the report, and you can transfer it to
your computer for viewing using Fugu or scp.

# References

## Kenya project - subtype consensuses for A, AE, C, & D

These subtype consensuses are based on sequences from the LANL database, and
then trimmed to the amplicon of interest for the Kenya project.  They are all
the same length, and use the same codon coordinates (BED files).  Additionally,
the amplicon was trimmed by one base at the beginning so that the CDS starts at
the first reference base instead of the second reference base.

The Kenya references are currently special because the pipeline knows how to
generate GenBank files for SnpEff from the fasta.

## AF033819_pol

refs/AF033819_pol.gb is downloaded from GenBank:

    http://www.ncbi.nlm.nih.gov/nuccore/AF033819.3?from=1385&to=4964&sat=4&sat_key=39317181&report=gbwithparts

It is a partial download of the HIV-1 genome AF033819.3 which completely covers
Pol and includes substantial upstream/downstream margins.

This was figured out by finding where the original reference sequence matched.

# Required tools

Various software packages and tools are used by the code in this repository.
Some of them are provided in and used from the tools/ directory, others are
expected to be installed on the system and available in your PATH.

Included in tools/

    • SnpEff and SnpSift
    • FastQC
    • various custom utilities

Must be in PATH

    • Oracle Java 1.7+
    • sickle (at least v1.33-20-g4b0dc85-mullins from https://github.com/MullinsLab/sickle)
    • bwa
    • samtools (at least version 1.1)
    • bgzip and tabix
    • LoFreq >= 2.0
    • freebayes (at least v0.9.20-22-gd9fc3ac from https://github.com/MullinsLab/freebayes)
    • vcflib (https://github.com/ekg/vcflib)
    • App::RecordStream - recs commands (from CPAN)
    • Perl 5.018 or newer
    • pysam 0.8.2.1 or newer

Perl modules installed from CPAN:

    • Module::Runtime
    • BioPerl
    • Bio::Cigar
    • App::RecordStream - recs commands
    • App::RecordStream::Bio

I suggest using
[cpanm](https://metacpan.org/pod/distribution/App-cpanminus/bin/cpanm) to
install dependencies easily and quickly.

It's possible, and recommended, to use something like
[perlbrew](http://perlbrew.pl) or
[plenv](https://github.com/tokuhirom/plenv#readme) to setup an isolated Perl
install for this pipeline.  The advantages of this are that it's decoupled from
the system Perl and doesn't require root privileges to install modules.  You'll
want to adjust `PERLBIN` in the `Makefile` if you do so and don't use the same
path.
