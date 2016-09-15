C=$(shell pwd)
SPACE :=
SPACE +=
COMMA := ,
NCPU=4
MININS=0
MAXINS=220

DATA=data
PREFIX=SN15
GENOME_FILE=$(DATA)/Parastagonospora_nodorum_SN15_scaffolds.fasta
ANNOTATION_FILE=$(DATA)/Parastagonospora_nodorum_SN15.gtf
ANNOTATION_BED=$(DATA)/Parastagonospora_nodorum_SN15.bed
READS1_PATTERN=-R1.fastq.gz
READS2_PATTERN=-R2.fastq.gz
READS1=$(shell ls $(DATA)/*$(READS1_PATTERN))
READS2=$(shell ls $(DATA)/*$(READS2_PATTERN))
#READS1=data/invitro-SN15-000-R1.fastq.gz
#READS2=data/invitro-SN15-000-R2.fastq.gz
SAMPLE_NAMES=$(subst $(READS1_PATTERN),,$(notdir $(READS1)))

#HISAT_DOCKER=sudo docker run --rm -v $$PWD:/data:z darcyabjones/hisat2
HISAT_DOCKER=
#SAMTOOLS_DOCKER=sudo docker run --rm -v $$PWD:/data:z darcyabjones/samtools
SAMTOOLS_DOCKER=
RSEQC_DOCKER=
MULTIQC_DOCKER=
HTSEQ_DOCKER=
#STRINGTIE_DOCKER=sudo docker run --rm -v $$PWD:/data:z darcyabjones/stringtie
STRINGTIE_DOCKER=

# 01 - build index

GENOME_INDEX_DIR=genome_index
GENOME_INDEX_FILES=$(foreach i, 1, $(GENOME_INDEX_DIR)/$(PREFIX).$i.ht2)

SPLICE_FILE=$(GENOME_INDEX_DIR)/splice_sites.txt
EXTRACT_SPLICE=$(HISAT_DOCKER) extract_splice_sites.py $(1) > $(SPLICE_FILE)
NOVEL_SPLICE_FILE=$(GENOME_INDEX_DIR)/novel_splice_sites.txt

# 02 - align to reference
# Complete ones have all reads (even unaligned, use this for QC), regular ones just have the mapped ones.

ALIGN_DIR=align

BAM_COMPLETE_EXTS=.complete.bam
BAM_COMPLETE_TARGETS=$(foreach e, $(BAM_COMPLETE_EXTS), $(addprefix $(ALIGN_DIR)/, $(addsuffix $(e), %)))
BAM_COMPLETE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(BAM_COMPLETE_TARGETS)))
BAI_COMPLETE_FILES=$(addsuffix .bai, $(BAM_COMPLETE_FILES))

BAM_FILES=$(addprefix $(ALIGN_DIR)/, $(addsuffix .bam, $(SAMPLE_NAMES)))
BAI_FILES=$(addsuffix .bai, $(BAM_FILES))

# 03 - QC

QC_DIR=bam_qc

BAM_STAT_FILES=$(addprefix $(QC_DIR)/, $(addsuffix .bam_stat, $(notdir $(BAM_FILES))))

CLIPPING_PROFILE_EXTS=-clipping_profile.r .clipping_profile.R1.pdf .clipping_profile.R2.pdf .clipping_profile.xls
CLIPPING_PROFILE_TARGETS=$(foreach e, $(CLIPPING_PROFILE_EXTS), $(addprefix $(QC_DIR)/, $(addsuffix $(e), %)))
CLIPPING_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(CLIPPING_PROFILE_TARGETS)))

DELETION_PROFILE_EXTS=.deletion_profile.r .deletion_profile.txt
DELETION_PROFILE_TARGETS=$(foreach e, $(DELETION_PROFILE_EXTS), $(addprefix $(QC_DIR)/, $(addsuffix $(e), %)))
DELETION_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(DELETION_PROFILE_TARGETS)))

INNER_DISTANCE_EXTS=.inner_distance_freq.txt .inner_distance_plot.pdf .inner_distance_plot.r .inner_distance_plot.r
INNER_DISTANCE_TARGETS=$(foreach e, $(INNER_DISTANCE_EXTS), $(addprefix $(QC_DIR)/, $(addsuffix $(e), %)))
INNER_DISTANCE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(INNER_DISTANCE_TARGETS)))

INSERTION_PROFILE_EXTS=.insertion_profile.r .insertion_profile.R1.pdf .insertion_profile.R2.pdf .insertion_profile.xls
INSERTION_PROFILE_TARGETS=$(foreach e, $(INSERTION_PROFILE_EXTS), $(addprefix $(QC_DIR)/, $(addsuffix $(e), %)))
INSERTION_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(INSERTION_PROFILE_TARGETS)))

JUNCTION_ANNOTATION_EXTS=.junction.bed .junction_plot.r .junction.xls
JUNCTION_ANNOTATION_TARGETS=$(foreach e, $(JUNCTION_ANNOTATION_EXTS), $(addprefix $(QC_DIR)/, $(addsuffix $(e), %)))
JUNCTION_ANNOTATION_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(JUNCTION_ANNOTATION_TARGETS)))

JUNCTION_SATURATION_EXTS=.junctionSaturation_plot.pdf .junctionSaturation_plot.r
JUNCTION_SATURATION_TARGETS=$(foreach e, $(JUNCTION_SATURATION_EXTS), $(addprefix $(QC_DIR)/, $(addsuffix $(e), %)))
JUNCTION_SATURATION_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(JUNCTION_SATURATION_TARGETS)))

MISMATCH_PROFILE_EXTS=.mismatch_profile.pdf .mismatch_profile.r .mismatch_profile.xls
MISMATCH_PROFILE_TARGETS=$(foreach e, $(MISMATCH_PROFILE_EXTS), $(addprefix $(QC_DIR)/, $(addsuffix $(e), %)))
MISMATCH_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(MISMATCH_PROFILE_TARGETS)))

READ_DUPLICATION_EXTS=.DupRate_plot.pdf .DupRate_plot.r .pos.DupRate.xls .seq.DupRate.xls
READ_DUPLICATION_TARGETS=$(foreach e, $(READ_DUPLICATION_EXTS), $(addprefix $(QC_DIR)/, $(addsuffix $(e), %)))
READ_DUPLICATION_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(READ_DUPLICATION_TARGETS)))

GENEBODY_COVERAGE_EXTS=.geneBodyCoverage.curves.pdf .geneBodyCoverage.r .geneBodyCoverage.txt
GENEBODY_COVERAGE_TARGETS=$(foreach e, $(GENEBODY_COVERAGE_EXTS), $(addsuffix $(e), $(QC_DIR)/%))
GENEBODY_COVERAGE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(GENEBODY_COVERAGE_TARGETS)))

QC_FILES=$(BAM_STAT_FILES) $(CLIPPING_PROFILE_FILES) \
  $(DELETION_PROFILE_FILES) $(INNER_DISTANCE_FILES) \
	$(INSERTION_PROFILE_FILES) $(JUNCTION_ANNOTATION_FILES) \
	$(JUNCTION_SATURATION_FILES) $(MISMATCH_PROFILE_FILES) \
	$(READ_DUPLICATION_FILES) $(GENEBODY_COVERAGE_FILE)

MULTIQC_FILES=$(QC_DIR)/multiqc_report.html

# 04 - Count reads

COUNT_DIR=counts
COUNT_FILES=$(addprefix $(COUNT_DIR)/, $(addsuffix -counts.tsv, $(SAMPLE_NAMES)))

# 05 - Calculate coverage

COVERAGE_DIR=coverage
COVERAGE_FILES=$(addprefix $(COVERAGE_DIR)/, $(addsuffix -coverage.bedgraph, $(SAMPLE_NAMES)))

# 06 - Assemble transcripts.

ASSEMBLE_DIR=assembly
ASSEMBLE_EXTS=.stringtie.gtf .stringtie.ctab
ASSEMBLE_TARGETS=$(foreach e, $(ASSEMBLE_EXTS), $(addprefix $(ASSEMBLE_DIR)/, $(addsuffix $(e), %)))
ASSEMBLE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(ASSEMBLE_TARGETS)))


# Do the actual work
phony:
	$(HISAT_DOCKER) ls -l data

cleanall: all clean
	@rm -rf $^

clean: $(BAM_COMPLETE_FILES) $(BAI_COMPLETE_FILES)
	@rm -rf $^

all: index align qc count assemble

index: $(GENOME_INDEX_FILES) $(SPLICE_FILE) $(NOVEL_SPLICE_FILE)
align: $(BAM_FILES) $(BAI_FILES)
qc: $(QC_FILES) $(MULTIQC_FILES)
count: $(COUNT_FILES)
coverage: $(COVERAGE_FILES)
assemble: $(ASSEMBLE_FILES)

# 01 - build index
$(GENOME_INDEX_FILES): $(GENOME_FILE)
	@mkdir -p $(dir $@)
	$(HISAT_DOCKER) hisat2-build -f $< $(GENOME_INDEX_DIR)/$(PREFIX)

$(SPLICE_FILE): $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	$(HISAT_DOCKER) extract_splice_sites.py $< > $@

$(NOVEL_SPLICE_FILE): $(SPLICE_FILE) $(GENOME_INDEX_FILES) $(READS1) $(READS2)
	@mkdir -p $(dir $@)
	$(HISAT_DOCKER) hisat2 \
	  --threads $(NCPU) \
	  --known-splicesite-infile $(SPLICE_FILE) \
		--minins $(MININS) \
		--maxins $(MAXINS) \
		--min-intronlen 20 \
		--max-intronlen 10000 \
		--rna-strandness FR \
		--no-mixed \
		--no-discordant \
		--no-temp-splicesite \
		--downstream-transcriptome-assembly \
		--novel-splicesite-outfile $@ \
		-x $(GENOME_INDEX_DIR)/$(PREFIX) \
		-1 $(subst $(SPACE),$(COMMA),$(READS1)) \
		-2 $(subst $(SPACE),$(COMMA),$(READS2)) \
		| $(SAMTOOLS_DOCKER) samtools view -bT $(GENOME_FILE) -o $(NOVEL_SPLICE_FILE).bam -

# 02 - align to reference
$(BAM_COMPLETE_TARGETS): $(DATA)/%$(READS1_PATTERN) $(DATA)/%$(READS2_PATTERN) $(SPLICE_FILE) $(NOVEL_SPLICE_FILE) $(GENOME_INDEX_FILES)
	@mkdir -p $(dir $@)
	$(HISAT_DOCKER) hisat2 \
	  --threads $(NCPU) \
	  --known-splicesite-infile $(SPLICE_FILE) \
		--novel-splicesite-infile $(NOVEL_SPLICE_FILE) \
		--minins $(MININS) \
		--maxins $(MAXINS) \
		--min-intronlen 20 \
		--max-intronlen 10000 \
		--downstream-transcriptome-assembly \
		--rna-strandness FR \
		--no-mixed \
		--no-discordant \
		-x $(GENOME_INDEX_DIR)/$(PREFIX) \
		-1 $(word 1, $^) \
		-2 $(word 2, $^) \
		-S $(basename $(basename $@)).complete.sam
	$(SAMTOOLS_DOCKER) samtools view \
	  	-uT $(GENOME_FILE) \
			$(basename $(basename $@)).complete.sam \
		| $(SAMTOOLS_DOCKER) samtools sort \
		  -O BAM \
			-@ $(NCPU) \
			-l 9 \
			-o $(basename $(basename $@)).complete.bam

	@rm $(basename $(basename $@)).complete.sam

$(ALIGN_DIR)/%.bam: $(ALIGN_DIR)/%.complete.bam
	$(SAMTOOLS_DOCKER) samtools view -bh -F 0x0004 -@ $(NCPU) $< > $@.tmp \
		&& mv $@.tmp $@

$(ALIGN_DIR)/%.bai: $(ALIGN_DIR)/%
	@mkdir -p $(dir $@)
	$(SAMTOOLS_DOCKER) samtools index $<

# 03 - QC

$(QC_DIR)/%.bam.bam_stat: $(ALIGN_DIR)/%.complete.bam $(ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) bam_stat.py -i $< > $@.tmp \
		&& mv $@.tmp $@

$(CLIPPING_PROFILE_TARGETS): $(ALIGN_DIR)/%.complete.bam $(ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) clipping_profile.py \
		-i $< \
		-s "PE" \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(DELETION_PROFILE_TARGETS): $(ALIGN_DIR)/%.complete.bam $(ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) deletion_profile.py \
		-i $< \
		-l 125 \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(INNER_DISTANCE_TARGETS): $(ALIGN_DIR)/%.complete.bam $(ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) inner_distance.py \
		-i $< \
		--refgene=$(ANNOTATION_BED) \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(INSERTION_PROFILE_TARGETS): $(ALIGN_DIR)/%.complete.bam $(ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) insertion_profile.py \
		-i $< \
		-s PE \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(JUNCTION_ANNOTATION_TARGETS): $(ALIGN_DIR)/%.complete.bam $(ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) junction_annotation.py \
		-i $< --refgene=$(ANNOTATION_BED) \
		--min-intron=5 \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(JUNCTION_SATURATION_TARGETS): $(ALIGN_DIR)/%.complete.bam $(ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) junction_saturation.py \
		-i $< \
		--refgene=$(ANNOTATION_BED) \
		--min-intron=5 \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(MISMATCH_PROFILE_TARGETS): $(ALIGN_DIR)/%.complete.bam $(ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) mismatch_profile.py \
		-i $< \
		-l 125 \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(READ_DUPLICATION_TARGETS): $(ALIGN_DIR)/%.complete.bam $(ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) read_duplication.py \
		-i $< \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(GENEBODY_COVERAGE_TARGETS): $(BAM_COMPLETE_FILES) $(BAI_COMPLETE_FILES)
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) geneBody_coverage.py \
		-i $(subst $(SPACE),$(COMMA),$(BAM_COMPLETE_FILES))
		--refgene=$(ANNOTATION_BED) \
		-o $(QC_DIR)/combined

$(QC_DIR)/multiqc_report.html: $(BAM_COMPLETE_FILES) $(QC_FILES)
	@mkdir -p $(dir $@)
	$(MULTIQC_DOCKER) multiqc -o $(dir $@) $(QC_DIR)

# 04 - count reads

# Note, because we sorted by position htseq-count uses a huge amount of RAM to count fragments.
# Inplanta counts used about 8 GB of RAM to count, invitro counts used 15 GB of RAM, if this is too much for your computer i suggest looking at a different program or temporarily sorting by name.
# Also by default, htseq-count sets a maximum ram useage of ~3GB, you need to edit the package following http://seqanswers.com/forums/showthread.php?p=197997 .

$(COUNT_DIR)/%-counts.tsv: $(ALIGN_DIR)/%.bam $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	$(HTSEQ_DOCKER) htseq-count --format=bam --order=pos --type=exon $(word 1, $^) $(word 2, $^) > $@.tmp \
	  && mv $@.tmp $@

# 05 - Coverage

COVERAGE_FILES=$(addprefix $(COVERAGE_DIR), $(addsuffix -coverage.bedgraph, $(SAMPLE_NAMES)))

$(COVERAGE_DIR)/%-coverage.bedgraph: $(ALIGN_DIR)/%.bam $(GENOME_FILE)
	@mkdir -p $(dir $@)
	bedtools genomecov -bga -split -trackline -ibam $(word 1, $^) -g $(word 2, $^) > $@.tmp \
	  && mv $@.tmp $@

# 06 - stringtie

$(ASSEMBLE_TARGETS): $(ALIGN_DIR)/%.bam $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	$(STRINGTIE_DOCKER) stringtie \
	  $(word 1, $^) \
	  -p $(NCPU) \
		-G $(ANNOTATION_FILE) \
		-e \
		-b $(dir $@)/$(notdir $(basename $(word 1, $^))) \
		-o $(dir $@)/$(notdir $(basename $(word 1, $^))).stringtie.gtf
