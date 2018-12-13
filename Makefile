SHELL:=/bin/bash
export SHELLOPTS:=errexit:pipefail
.DELETE_ON_ERROR:
.SECONDARY:

TOOLS  = ./tools/
LOFREQ = lofreq
SNPEFF = $(TOOLS)/snpEff
FastQC = $(TOOLS)/FastQC/fastqc
BWAOPT = -A 1 -B 1 -O 1 -E 1
FREEBAYES = freebayes

MIN_TAG_FAMILY_SIZE = 3

TRIMMOMATIC          = java -jar $(TOOLS)/Trimmomatic/trimmomatic.jar PE -phred33 -threads 1
TRIM_ILLUMINACLIP    = $(TRIM_SEED_MISMATCHES):$(TRIM_PALINDROME):$(TRIM_SIMPLE)
TRIM_ADAPTERS        = $(TOOLS)/Trimmomatic/adapters/NexteraPE-PE.fa
TRIM_SEED_MISMATCHES = 2
TRIM_SIMPLE          = 10
TRIM_PALINDROME      = 30

TRIM5PRIMER = $(TOOLS)/Trim5Primer

PERLBIN := /opt/perl/5.20/bin
export PATH := $(PERLBIN):/usr/local/bin:$(PATH)

# Any (non-phony) target using this should declare a dependency on amplicon-ref-by-run.csv
REFS = $(shell recs fromcsv --header amplicon-ref-by-run.csv | recs eval '{{reference}}' | sort -u)

# $(call if-tagged-run,...)
# $1 = Run ID
# $2 = string if run has template tags (true condition)
# $3 = string if run lacks template tags (false condition)
# $4 = substitution pattern
# $5 = string which is returned after replacing $4 with $2 or $3
#
# This function is most likely useful with secondary expansion and called like:
#	$$(call if-tagged-run,$$*,A,B,{},cache/{}/%.bam)
#
if-tagged-run = $(if $(shell $(TOOLS)/run-is-tagged $1),$(subst $4,$2,$5),$(subst $4,$3,$5))

runrefname = $(shell $(TOOLS)/amplicon-ref-for-run $(1))
runref     = refs/$(call runrefname,$(1)).fasta
runbed     = refs/$(call runrefname,$(1)).bed

FASTQGZ = $(wildcard reads/*.fastq.gz)
RUNS    = $(filter-out $(EXCLUDE_RUNS),$(patsubst reads/%_R1.fastq.gz,%,$(filter %_R1.fastq.gz,$(FASTQGZ))))
EXCLUDE_RUNS =

#
# Generic targets
#

# make directories used by order-only prerequisites, e.g.
# 	cache/somedir/output: input | mkcache-somedir
#
mkcache-%:
	test -d cache/$* || mkdir -p cache/$*

mkdir-%:
	test -d $* || mkdir -p $*

%.fasta.fai: %.fasta
	samtools faidx $<

%.fasta.amb %.fasta.ann %.fasta.bwt %.fasta.pac %.fasta.sa: %.fasta
	bwa index $<

%.bam: %.sam
	samtools view -b $< > $@

%.bam.bai: %.bam
	samtools index $<

%.vcf.gz %.vcf.gz.tbi: %.vcf
	bgzip -c $< > $@
	tabix -p vcf $@

define dofastqc
	$(FastQC) $<
	mv $*{_,.}fastqc.html
	rm -v $*_fastqc.zip
endef

%.fastqc.html: %.fastq
	$(dofastqc)
%.fastqc.html: %.fastq.gz
	$(dofastqc)

#
# Analysis-specific targets
#

.SECONDEXPANSION: cache/ref/%.fasta cache/ref/%.bed cache/ref/%.primers.fasta
cache/ref/%.fasta: $$(call runref,$$*) | mkcache-ref
	cp $< $@
cache/ref/%.bed: $$(call runbed,$$*) cache/ref/%.fasta | mkcache-ref
	cp $< $@
	$(TOOLS)/check-bed $^

cache/ref/%.primers.fasta: refs/$$(call runrefname,$$*).primers.fasta | mkcache-ref
	recs fromfasta --oneline $< \
		| recs assert -v -d 'Primer does not contain ambiguity codes' '{{seq}} =~ /^[ATCG]+$$/' >/dev/null
	cp $< $@

cache/ref/%.clip.fasta: cache/ref/%.primers.fasta $(TRIM_ADAPTERS) $(TOOLS)/rc-primers.xform | mkcache-ref
	cp $(TRIM_ADAPTERS) $@
	# Adapter clipping can only use the RC of primers since it looks for 3'
	# contamination.  The RC's remain limited to a single read direction
	# (opposite what the original primer was, of course).
	recs fromfasta --oneline $< \
		| recs grep '{{id}} =~ m{/[12]$$}' \
		| recs xform -E $(TOOLS)/rc-primers.xform \
		| recs tofasta >> $@

# Adapter clipping
$(foreach out,R1 R2 R1_orphans R2_orphans,cache/clipped/%_$(out).fastq.gz): reads/%_R1.fastq.gz reads/%_R2.fastq.gz cache/ref/%.clip.fasta | mkcache-clipped
	$(TRIMMOMATIC) -trimlog cache/clipped/$*.log \
		reads/$*_R[1-2].fastq.gz \
		cache/clipped/$*_{R1{,_orphans},R2{,_orphans}}.fastq.gz \
		ILLUMINACLIP:cache/ref/$*.clip.fasta:$(TRIM_ILLUMINACLIP)

# Tag extraction
cache/tagged/%_R1.fastq.gz cache/tagged/%_R2.fastq.gz cache/tagged/%.tags cache/tagged/%.badtags: $(TOOLS)/rt-tag-to-header cache/clipped/%_R1.fastq.gz cache/clipped/%_R2.fastq.gz | mkcache-tagged
	$(TOOLS)/rt-tag-to-header \
		--length 8 \
		--output-tags cache/tagged/$*.tags \
		--output-bad-tags cache/tagged/$*.badtags \
		--homopolymer-threshold 9 \
		cache/{clipped,tagged}/$*_R{1,2}.fastq.gz

# These summary stats and plots on tags are used to examine/QC the raw observed
# tag list from the FastQ files as well as the tags retained after mapping in
# the cache/tag-groups/ targets.
%.tagcounts: %.tags
	sort $< | uniq -c | perl -anE 'say join "\t", @F[1,0]' > $@

%.tagstats: %.tagcounts
	# Similar to ConsensusMaker.py:tagStats(), but without the tag family size
	# cutoff applying to cigar subgroups.
	recs fromcsv -d $$'\t' -k tag,size $< \
		| recs grep '{{size}} >= $(MIN_TAG_FAMILY_SIZE)' \
		| recs collate -k size -a count \
		| recs xform '{{count}} *= {{size}}' \
		| recs collate -a total=sum,count -a sizes=records \
		| recs xform '$$_->{count} /= {{total}}, push_output($$_) for @{ {{sizes}} }' \
		| recs tocsv -d $$'\t' -k size,count --noheader \
		> $@

%.template-counts: %.tagcounts
	recs fromcsv -d $$'\t' -k tag,size $< \
		| recs grep '{{size}} >= $(MIN_TAG_FAMILY_SIZE)' \
		| recs collate -k size -a tag_count=count \
		| recs sort -k size=numeric \
		| recs tocsv -d $$'\t' -k size,tag_count \
		> $@

%.template-counts.png: $(TOOLS)/plot-templates-by-family-size %.template-counts
	$^ $@ $(MIN_TAG_FAMILY_SIZE)

%.family-sizes.png: $(TOOLS)/plot-family-size %.tagstats
	$^ $@ $(MIN_TAG_FAMILY_SIZE)

# 5' primer clipping
.SECONDEXPANSION: cache/clipped5p/%_R1.fastq.gz cache/clipped5p/%_R2.fastq.gz
cache/clipped5p/%_R1.fastq.gz cache/clipped5p/%_R2.fastq.gz: cache/ref/%.primers.fasta $$(call if-tagged-run,$$*,tagged,clipped,{},cache/{}/$$*_R1.fastq.gz cache/{}/$$*_R2.fastq.gz) | mkcache-clipped5p
	$(TRIM5PRIMER)/trim5primer.py --primers $^ cache/clipped5p/$*_R{1,2}.fastq.gz

# Quality control
cache/qc/%_R1.fastq cache/qc/%_R2.fastq cache/qc/%_orphans.fastq: cache/clipped5p/%_R1.fastq.gz cache/clipped5p/%_R2.fastq.gz | mkcache-qc
	# "sanger" just means Phred+33 which Illumina (CASAVA 1.8) uses.
	sickle pe \
		--qual-type=sanger \
		--qual-threshold=30 \
		--length-threshold=75 \
		--window=9 \
		--truncate-n \
		-f <(gunzip -c $(word 1,$^)) \
		-r <(gunzip -c $(word 2,$^)) \
		-o cache/qc/$*_R1.fastq \
		-p cache/qc/$*_R2.fastq \
		-s cache/qc/$*_orphans.fastq \
		| tee cache/qc/$*.log

cache/qc/%.stats: $(TOOLS)/sickle-summary reads/%_R1.fastq.gz reads/%_R2.fastq.gz cache/qc/%_R1.fastq cache/qc/%_R2.fastq cache/qc/%_orphans.fastq | mkcache-qc
	$^ > $@

qc: $(patsubst %,cache/qc/%_R1.fastq,$(RUNS))

# Mapping, against a run-specific reference
cache/mapped/%.sam: cache/qc/%_R1.fastq cache/qc/%_R2.fastq cache/ref/%.fasta cache/ref/%.fasta.bwt | mkcache-mapped
	bwa mem -t 1 $(BWAOPT) \
		-R '@RG\tID:1\tPL:Illumina\tPU:unknown\tLB:Nextera\tSM:$*' \
		-M cache/ref/$*.fasta \
		cache/qc/$*_R{1,2}.fastq > $@

mapped: $(patsubst %,cache/mapped/%.sam,$(RUNS))

# 
# Tagged reads -> Consensus sequences
#
cache/tag-groups/%.tags cache/tag-groups/%.filtered.tags: cache/mapped/%.sam $(TOOLS)/sam-filter-to-overlapping-pairs | mkcache-tag-groups
	# Forward reads of fully mapped pairs, no secondary or supplementary
	# alignments or bad QC.  Extract tags from reads, filter to tags
	# observed for MIN_TAG_FAMILY_SIZE or more read pairs.
	samtools view -f 0x43 -F 0xC0C $< \
		| $(TOOLS)/sam-filter-to-overlapping-pairs \
		| perl -anE 'say $$1 if $$F[0] =~ /\|([ATCG]+)$$/' \
		| tee cache/tag-groups/$*.tags \
		| recs fromcsv -k tag \
		| recs collate -k tag -a count \
		| recs grep '{{count}} >= $(MIN_TAG_FAMILY_SIZE)' \
		| recs eval '{{tag}}' \
		> cache/tag-groups/$*.filtered.tags

cache/tag-groups/%.bam: $(TOOLS)/tags-to-read-groups cache/tag-groups/%.filtered.tags cache/mapped/%.sam | mkcache-tag-groups
	$^ /dev/stdout | samtools sort -O bam -T $@ -o $@ /dev/stdin

# Call a separate consensus for each tag group
cache/consensus/%.fq: cache/tag-groups/%.bam cache/tag-groups/%.bam.bai $(TOOLS)/simple-consensus-per-read-group | mkcache-consensus
	$(TOOLS)/simple-consensus-per-read-group --threshold 0.7 $< > $@

# Map consensus seqs to reference again, now unpaired
cache/consensus/%.sam: cache/ref/%.fasta cache/consensus/%.fq | mkcache-consensus
	bwa mem -t 1 $(BWAOPT) -R '@RG\tID:1\tPL:Illumina\tPU:unknown\tLB:Nextera\tSM:$*' -M $^ > $@

# Filter to only reads intersecting known DRM codons
.SECONDEXPANSION: cache/drm-codons/%.bam
cache/drm-codons/%.bam: $$(call if-tagged-run,$$*,cache/consensus/$$*.sam,cache/mapped/$$*.sam,x,x) cache/ref/%.fasta cache/ref/%.bed | mkcache-drm-codons
	samtools view -L cache/ref/$*.bed -uT cache/ref/$*.fasta $< \
		| samtools sort -O bam -T $*.drm-codons \
		> $@

drm-codons: $(patsubst %,cache/drm-codons/%.bam,$(RUNS))

cache/lofreq/%.vcf: cache/drm-codons/%.bam cache/ref/%.fasta cache/ref/%.bed | mkcache-lofreq
	rm -f $@.tmp
	time $(LOFREQ) call -a 1 --no-default-filter \
		--bed cache/ref/$*.bed \
		--ref cache/ref/$*.fasta \
		--out $@.tmp \
		$<
	mv -f $@.tmp $@

lofreq: $(patsubst %,cache/lofreq/%.vcf,$(RUNS))

cache/haplotypes/%.vcf: cache/lofreq/%.vcf.gz cache/drm-codons/%.bam cache/drm-codons/%.bam.bai cache/ref/%.fasta | mkcache-haplotypes
	time $(FREEBAYES) \
		-f cache/ref/$*.fasta \
		--haplotype-basis-alleles cache/lofreq/$*.vcf.gz \
		--report-all-haplotype-alleles \
		--haplotype-length 1 \
		--pooled-continuous \
		--min-alternate-fraction 0 \
		--no-population-priors --hwe-priors-off --binomial-obs-priors-off --allele-balance-priors-off \
		cache/drm-codons/$*.bam \
		> $@

haplotypes: $(patsubst %,cache/haplotypes/%.vcf,$(RUNS))

# SnpEff - gather metadata about SNPs
snpEffDatabases: $(patsubst %,$(SNPEFF)/data/%/snpEffectPredictor.bin,$(REFS))

snpeff-config-%: snpEff.config
	# Configures SnpEff to know about the data/%/ directory for the reference "genome"
	grep -q "$*.genome" $< \
		|| echo "$*.genome : HIV-1" >> $<

$(SNPEFF)/data/%/genes.gbk: refs/%.gb
	[ -d $(SNPEFF)/data/$*/ ] || mkdir -p $(SNPEFF)/data/$*/
	cp $< $@

$(SNPEFF)/data/kenya_%/genes.gbk: refs/kenya_%.fasta
	[ -d $(SNPEFF)/data/kenya_$*/ ] || mkdir -p $(SNPEFF)/data/kenya_$*/
	$(TOOLS)/fake-genbank-for-pol-ref < $< > $@

$(SNPEFF)/data/PEPFAR_RT%_8E5/genes.gbk: refs/PEPFAR_RT%_8E5.fasta
	[ -d $(SNPEFF)/data/PEPFAR_RT$*_8E5/ ] || mkdir -p $(SNPEFF)/data/PEPFAR_RT$*_8E5/
	$(TOOLS)/fake-genbank-for-pol-ref < $< > $@

$(SNPEFF)/data/%/snpEffectPredictor.bin: $(SNPEFF)/data/%/genes.gbk | snpeff-config-%
	java -jar $(SNPEFF)/snpEff.jar build -v -genbank $*

cache/snpeff/%.vcf: cache/haplotypes/%.vcf snpEffDatabases | mkcache-snpeff
	vcfbreakmulti $< | java -jar $(SNPEFF)/snpEff.jar eff \
		-v -no-upstream -no-downstream -noStats -formatEff \
		$(call runrefname,$*) \
		> $@

cache/snp-codons/%.vcf: cache/snpeff/%.vcf cache/ref/%.bed | mkcache-snp-codons
	vcfannotate --bed cache/ref/$*.bed --key CODON_DRM_NAME $< > $@

snp-codons: $(patsubst %,cache/snp-codons/%.vcf,$(RUNS))

results/%.tsv: cache/snp-codons/%.vcf | mkdir-results
	vcf2tsv $< \
		| recs-fromcsv -d $$'\t' --header \
		| recs-xform '{{AF}}  = sprintf "%0.6f", {{AO}}  / {{DP}};' \
		| recs-xform '{{SAF}} = sprintf "%0.6f", {{SAF}} / {{DP}};' \
		| recs-xform '{{SAR}} = sprintf "%0.6f", {{SAR}} / {{DP}};' \
		| recs-xform '({{REF_CODON}}, {{ALT_CODON}}) = split /\//, (split /\|/, {{EFF}})[2];' \
		| recs-xform '({{REF_AA}},    {{ALT_AA}})    = (split /\|/, {{EFF}})[3] =~ m{^p\.([A-Z]+|[*])\d+([A-Z]+|[*])/}i;' \
		| recs-tocsv -d $$'\t' \
			-k CHROM,POS,ID,REF,ALT,QUAL,FILTER,CODON_DRM_NAME,EFF,REF_CODON,ALT_CODON,REF_AA,ALT_AA \
			-k `$(TOOLS)/vcfinfofields -d, $<` \
		> $@

summary/%.tsv: results/%.tsv | mkdir-summary
	recs fromcsv -d $$'\t' --header $< \
		| recs tocsv -d $$'\t' \
			-k 'CHROM,POS,REF,ALT,REF_CODON,ALT_CODON,REF_AA,ALT_AA,CODON_DRM_NAME,AF,DP,MQM,SAF,SAR' \
	> $@

all-results: $(patsubst %,results/%.tsv,$(RUNS))
all-summary: $(patsubst %,summary/%.tsv,$(RUNS))

tag_plots := $(foreach step,tagged tag-groups,cache/$(step)/%.template-counts.png cache/$(step)/%.family-sizes.png)
.SECONDEXPANSION: %
$(RUNS): %: summary/%.tsv $$(call if-tagged-run,$$*,$(tag_plots),,x,x)

.PHONY: clean
clean:
	rm -rfv cache/* results/*.tsv summary/*.tsv $(FASTQC_HTML)

versions:
	@printf '\n\n### PATH = %s\n' $$PATH
	@printf '\n\n### %s\n\n' `which samtools`
	samtools --version
	@printf '\n\n### %s\n\n' `which bwa`
	bwa || true
	@printf '\n\n### %s\n\n' `which perl`
	perl -V
	@printf '\n\n### %s\n\n' `which recs`
	recs --version
	@printf '\n\n### %s\n\n' `which python`
	python --version
	@printf '\n\n### %s\n\n' `which R`
	R --version
	@printf '\n\n### %s\n\n' `which java`
	java -version
