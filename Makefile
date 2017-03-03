C=$(shell pwd)
SPACE :=
SPACE +=
COMMA := ,
NCPU=24

define NEWLINE


endef

DATA=$(C)/data
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
#STAR_DOCKER=sudo docker run --rm -v $PWD:/data:z darcyabjones/star
STAR_DOCKER=
#SAMTOOLS_DOCKER=sudo docker run --rm -v $$PWD:/data:z darcyabjones/samtools
SAMTOOLS_DOCKER=
RSEQC_DOCKER=
MULTIQC_DOCKER=
HTSEQ_DOCKER=
#STRINGTIE_DOCKER=sudo docker run --rm -v $$PWD:/data:z darcyabjones/stringtie
STRINGTIE_DOCKER=
SUBREAD_DOCKER=

# 01 - build index

GENOME_INDEX_DIR=genome_index

HISAT_GENOME_INDEX_DIR=$(GENOME_INDEX_DIR)/hisat2
HISAT_GENOME_INDEX_FILES=$(foreach i, 1, $(HISAT_GENOME_INDEX_DIR)/$(PREFIX).$i.ht2)
HISAT_SPLICE_FILE=$(HISAT_GENOME_INDEX_DIR)/splice_sites.txt
HISAT_NOVEL_SPLICE_FILE=$(HISAT_GENOME_INDEX_DIR)/novel_splice_sites.txt

STAR_GENOME_INDEX_DIR=$(GENOME_INDEX_DIR)/star
#STAR_GENOME_INDEX_FILES=$(addprefix $(STAR_GENOME_INDEX_DIR)/, \
#  chrLength.txt \
#  chrNameLength.txt \
#  chrName.txt \
#  chrStart.txt \
#  exonGeTrInfo.tab \
#  exonInfo.tab \
#  geneInfo.tab \
#  Genome \
#  genomeParameters.txt \
#  SA \
#  SAindex \
#  sjdbInfo.txt \
#  sjdbList.fromGTF.out.tab \
#  sjdbList.out.tab \
#  transcriptInfo.tab \
#)
STAR_GENOME_INDEX_FILES=$(addprefix $(STAR_GENOME_INDEX_DIR)/, SA)

STAR_NOVEL_SPLICE_FILES=$(addprefix $(STAR_GENOME_INDEX_DIR)/, $(addsuffix .SJ.out.tab, $(SAMPLE_NAMES)))

# 02 - align to reference
# Complete ones have all reads (even unaligned, use this for QC), regular ones just have the mapped ones.

ALIGN_DIR=align
HISAT_ALIGN_DIR=$(ALIGN_DIR)/hisat2
HISAT_BAM_COMPLETE_TARGETS=$(addprefix $(HISAT_ALIGN_DIR)/, $(addsuffix .complete.bam, %))
HISAT_BAM_COMPLETE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_BAM_COMPLETE_TARGETS)))
HISAT_BAI_COMPLETE_FILES=$(addsuffix .bai, $(HISAT_BAM_COMPLETE_FILES))

HISAT_BAM_FILES=$(addprefix $(HISAT_ALIGN_DIR)/, $(addsuffix .bam, $(SAMPLE_NAMES)))
HISAT_BAI_FILES=$(addsuffix .bai, $(HISAT_BAM_FILES))

STAR_ALIGN_DIR=$(ALIGN_DIR)/star
STAR_BAM_FILES=$(addprefix $(STAR_ALIGN_DIR)/, $(addsuffix .bam, $(SAMPLE_NAMES)))
STAR_BAI_FILES=$(addsuffix .bai, $(STAR_BAM_FILES))

# 03 - QC

QC_DIR=bam_qc
HISAT_QC_DIR=$(QC_DIR)/hisat2
STAR_QC_DIR=$(QC_DIR)/star

HISAT_BAM_STAT_FILES=$(addprefix $(HISAT_QC_DIR)/, $(addsuffix .bam_stat, $(notdir $(HISAT_BAM_FILES))))
STAR_BAM_STAT_FILES=$(addprefix $(STAR_QC_DIR)/, $(addsuffix .bam_stat, $(notdir $(STAR_BAM_FILES))))

CLIPPING_PROFILE_EXTS=.clipping_profile.r .clipping_profile.R1.pdf .clipping_profile.R2.pdf .clipping_profile.xls
HISAT_CLIPPING_PROFILE_TARGETS=$(foreach e, $(CLIPPING_PROFILE_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_CLIPPING_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_CLIPPING_PROFILE_TARGETS)))
STAR_CLIPPING_PROFILE_TARGETS=$(foreach e, $(CLIPPING_PROFILE_EXTS), $(addprefix $(STAR_QC_DIR)/, $(addsuffix $(e), %)))
STAR_CLIPPING_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(STAR_CLIPPING_PROFILE_TARGETS)))

DELETION_PROFILE_EXTS=.deletion_profile.r .deletion_profile.txt
HISAT_DELETION_PROFILE_TARGETS=$(foreach e, $(DELETION_PROFILE_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_DELETION_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_DELETION_PROFILE_TARGETS)))
STAR_DELETION_PROFILE_TARGETS=$(foreach e, $(DELETION_PROFILE_EXTS), $(addprefix $(STAR_QC_DIR)/, $(addsuffix $(e), %)))
STAR_DELETION_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(STAR_DELETION_PROFILE_TARGETS)))

INNER_DISTANCE_EXTS=.inner_distance_freq.txt .inner_distance_plot.pdf .inner_distance_plot.r .inner_distance_plot.r
HISAT_INNER_DISTANCE_TARGETS=$(foreach e, $(INNER_DISTANCE_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_INNER_DISTANCE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_INNER_DISTANCE_TARGETS)))
STAR_INNER_DISTANCE_TARGETS=$(foreach e, $(INNER_DISTANCE_EXTS), $(addprefix $(STAR_QC_DIR)/, $(addsuffix $(e), %)))
STAR_INNER_DISTANCE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(STAR_INNER_DISTANCE_TARGETS)))

INSERTION_PROFILE_EXTS=.insertion_profile.r .insertion_profile.R1.pdf .insertion_profile.R2.pdf .insertion_profile.xls
HISAT_INSERTION_PROFILE_TARGETS=$(foreach e, $(INSERTION_PROFILE_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_INSERTION_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_INSERTION_PROFILE_TARGETS)))
STAR_INSERTION_PROFILE_TARGETS=$(foreach e, $(INSERTION_PROFILE_EXTS), $(addprefix $(STAR_QC_DIR)/, $(addsuffix $(e), %)))
STAR_INSERTION_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(STAR_INSERTION_PROFILE_TARGETS)))

JUNCTION_ANNOTATION_EXTS=.junction.bed .junction_plot.r .junction.xls
HISAT_JUNCTION_ANNOTATION_TARGETS=$(foreach e, $(JUNCTION_ANNOTATION_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_JUNCTION_ANNOTATION_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_JUNCTION_ANNOTATION_TARGETS)))
STAR_JUNCTION_ANNOTATION_TARGETS=$(foreach e, $(JUNCTION_ANNOTATION_EXTS), $(addprefix $(STAR_QC_DIR)/, $(addsuffix $(e), %)))
STAR_JUNCTION_ANNOTATION_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(STAR_JUNCTION_ANNOTATION_TARGETS)))

JUNCTION_SATURATION_EXTS=.junctionSaturation_plot.pdf .junctionSaturation_plot.r
HISAT_JUNCTION_SATURATION_TARGETS=$(foreach e, $(JUNCTION_SATURATION_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_JUNCTION_SATURATION_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_JUNCTION_SATURATION_TARGETS)))
STAR_JUNCTION_SATURATION_TARGETS=$(foreach e, $(JUNCTION_SATURATION_EXTS), $(addprefix $(STAR_QC_DIR)/, $(addsuffix $(e), %)))
STAR_JUNCTION_SATURATION_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(STAR_JUNCTION_SATURATION_TARGETS)))

MISMATCH_PROFILE_EXTS=.mismatch_profile.pdf .mismatch_profile.r .mismatch_profile.xls
HISAT_MISMATCH_PROFILE_TARGETS=$(foreach e, $(MISMATCH_PROFILE_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_MISMATCH_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_MISMATCH_PROFILE_TARGETS)))
STAR_MISMATCH_PROFILE_TARGETS=$(foreach e, $(MISMATCH_PROFILE_EXTS), $(addprefix $(STAR_QC_DIR)/, $(addsuffix $(e), %)))
STAR_MISMATCH_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(STAR_MISMATCH_PROFILE_TARGETS)))

READ_DUPLICATION_EXTS=.DupRate_plot.pdf .DupRate_plot.r .pos.DupRate.xls .seq.DupRate.xls
HISAT_READ_DUPLICATION_TARGETS=$(foreach e, $(READ_DUPLICATION_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_READ_DUPLICATION_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_READ_DUPLICATION_TARGETS)))
STAR_READ_DUPLICATION_TARGETS=$(foreach e, $(READ_DUPLICATION_EXTS), $(addprefix $(STAR_QC_DIR)/, $(addsuffix $(e), %)))
STAR_READ_DUPLICATION_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(STAR_READ_DUPLICATION_TARGETS)))

GENEBODY_COVERAGE_EXTS=.geneBodyCoverage.curves.pdf .geneBodyCoverage.r .geneBodyCoverage.txt
HISAT_GENEBODY_COVERAGE_TARGETS=$(foreach e, $(GENEBODY_COVERAGE_EXTS), $(addsuffix $(e), $(HISAT_QC_DIR)/%))
HISAT_GENEBODY_COVERAGE_FILES=$(foreach s, combined, $(subst %,$(s), $(HISAT_GENEBODY_COVERAGE_TARGETS)))
STAR_GENEBODY_COVERAGE_TARGETS=$(foreach e, $(GENEBODY_COVERAGE_EXTS), $(addsuffix $(e), $(STAR_QC_DIR)/%))
STAR_GENEBODY_COVERAGE_FILES=$(foreach s, combined, $(subst %,$(s), $(STAR_GENEBODY_COVERAGE_TARGETS)))

HISAT_MULTIQC_FILES=$(HISAT_QC_DIR)/multiqc_report.html
STAR_MULTIQC_FILES=$(STAR_QC_DIR)/multiqc_report.html

HISAT_QC_FILES=$(HISAT_BAM_STAT_FILES) $(HISAT_CLIPPING_PROFILE_FILES) \
  $(HISAT_DELETION_PROFILE_FILES) $(HISAT_INNER_DISTANCE_FILES) \
	$(HISAT_INSERTION_PROFILE_FILES) $(HISAT_JUNCTION_ANNOTATION_FILES) \
	$(HISAT_JUNCTION_SATURATION_FILES) $(HISAT_MISMATCH_PROFILE_FILES) \
	$(HISAT_READ_DUPLICATION_FILES) $(HISAT_GENEBODY_COVERAGE_FILES)

STAR_QC_FILES=$(STAR_BAM_STAT_FILES) $(STAR_CLIPPING_PROFILE_FILES) \
	$(STAR_DELETION_PROFILE_FILES) $(STAR_INNER_DISTANCE_FILES) \
	$(STAR_INSERTION_PROFILE_FILES) $(STAR_JUNCTION_ANNOTATION_FILES) \
	$(STAR_JUNCTION_SATURATION_FILES) $(STAR_MISMATCH_PROFILE_FILES) \
	$(STAR_READ_DUPLICATION_FILES) $(STAR_GENEBODY_COVERAGE_FILES)

# 04 - Count reads

COUNT_DIR=counts

HISAT_COUNT_DIR=$(COUNT_DIR)/hisat2
HISAT_COUNT_FILES=$(addprefix $(HISAT_COUNT_DIR)/, feature_counts.tsv)

STAR_COUNT_DIR=$(COUNT_DIR)/star
STAR_COUNT_FILES=$(addprefix $(STAR_COUNT_DIR)/, feature_counts.tsv)

# 05 - Calculate coverage

COVERAGE_DIR=coverage
COVERAGE_EXTS=-coverage.bedgraph -coverage-fwd.bedgraph -coverage-rev.bedgraph

HISAT_COVERAGE_DIR=$(COVERAGE_DIR)/hisat2
HISAT_COVERAGE_FILES=$(foreach e, $(COVERAGE_EXTS), $(addprefix $(HISAT_COVERAGE_DIR)/, $(addsuffix $(e), $(SAMPLE_NAMES))))

STAR_COVERAGE_DIR=$(COVERAGE_DIR)/star
STAR_COVERAGE_FILES=$(foreach e, $(COVERAGE_EXTS), $(addprefix $(STAR_COVERAGE_DIR)/, $(addsuffix $(e), $(SAMPLE_NAMES))))

# 06a - Assemble transcripts. stringtie a

CUFFLINKS_ASS_DIR=cufflinks_assemble
CUFFLINKS_ASS_FILES=transcripts.gtf # isoforms.fpkm_tracking genes.fpkm_tracking

HISAT_CUFFLINKS_ASS_DIR=$(CUFFLINKS_ASS_DIR)/hisat2
HISAT_CUFFLINKS_ASS_TARGETS=$(foreach e, $(CUFFLINKS_ASS_FILES), $(addprefix $(HISAT_CUFFLINKS_ASS_DIR)/%/, $(e)))
HISAT_CUFFLINKS_ASS_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_CUFFLINKS_ASS_TARGETS)))

STAR_CUFFLINKS_ASS_DIR=$(CUFFLINKS_ASS_DIR)/star
STAR_CUFFLINKS_ASS_TARGETS=$(foreach e, $(CUFFLINKS_ASS_FILES), $(addprefix $(STAR_CUFFLINKS_ASS_DIR)/%/, $(e)))
STAR_CUFFLINKS_ASS_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(STAR_CUFFLINKS_ASS_TARGETS)))

# 06b - merge gtf files

HISAT_CUFFLINKS_ASS_MERGED=$(HISAT_CUFFLINKS_ASS_DIR)/merged.gtf
STAR_CUFFLINKS_ASS_MERGED=$(STAR_CUFFLINKS_ASS_DIR)/merged.gtf

# 07- count cufflinks

CUFF_COUNT_DIR=cufflinks_counts

HISAT_CUFF_COUNT_DIR=$(CUFF_COUNT_DIR)/hisat2
HISAT_CUFF_COUNT_FILES=$(addprefix $(HISAT_CUFF_COUNT_DIR)/, feature_counts.tsv)

STAR_CUFF_COUNT_DIR=$(CUFF_COUNT_DIR)/star
STAR_CUFF_COUNT_FILES=$(addprefix $(STAR_CUFF_COUNT_DIR)/, feature_counts.tsv)

# 08 cuff quant

ORIG_CUFF_QUANT_DIR=cufflinks_quant_orig
CUFF_QUANT_FILE=abundances.cxb
ORIG_STAR_CUFF_QUANT_DIR=$(ORIG_CUFF_QUANT_DIR)/star

ORIG_STAR_CUFF_QUANT_TARGETS=$(foreach e,$(CUFF_QUANT_FILE),$(addprefix $(ORIG_STAR_CUFF_QUANT_DIR)/%/,$(e)))
ORIG_STAR_CUFF_QUANT_FILES=$(foreach s,$(SAMPLE_NAMES),$(subst %,$(s),$(ORIG_STAR_CUFF_QUANT_TARGETS)))

NEW_CUFF_QUANT_DIR=cufflinks_quant_new
NEW_STAR_CUFF_QUANT_DIR=$(NEW_CUFF_QUANT_DIR)/star

NEW_STAR_CUFF_QUANT_TARGETS=$(foreach e,$(CUFF_QUANT_FILE),$(addprefix $(NEW_STAR_CUFF_QUANT_DIR)/%/,$(e)))
NEW_STAR_CUFF_QUANT_FILES=$(foreach s,$(SAMPLE_NAMES),$(subst %,$(s),$(NEW_STAR_CUFF_QUANT_TARGETS)))


# 09 cuffnorm

CUFF_NORM_FILE=genes.fpkm_table
ORIG_CUFF_NORM_DIR=cufflinks_norm_orig

ORIG_STAR_CUFF_NORM_DIR=$(ORIG_CUFF_NORM_DIR)/star
ORIG_STAR_CUFF_NORM_FILES=$(foreach e, $(CUFF_NORM_FILE), $(addprefix $(ORIG_STAR_CUFF_NORM_DIR)/, $(e)))


NEW_CUFF_NORM_DIR=cufflinks_norm_new
NEW_STAR_CUFF_NORM_DIR=$(NEW_CUFF_NORM_DIR)/star
NEW_STAR_CUFF_NORM_FILES=$(foreach e, $(CUFF_NORM_FILE), $(addprefix $(NEW_STAR_CUFF_NORM_DIR)/, $(e)))

phony:
	@echo $(subst $(SPACE),$(COMMA),$(NEW_STAR_CUFF_QUANT_FILES))
	@echo $(COMMA)

cleanall: all clean
	@rm -rf $^

clean: $(HISAT_BAM_COMPLETE_FILES) $(HISAT_BAI_COMPLETE_FILES)
	@rm -rf $^

all: index align qc count coverage

hisat_index: $(HISAT_GENOME_INDEX_FILES) $(HISAT_SPLICE_FILE) $(HISAT_NOVEL_SPLICE_FILE)
star_index: $(STAR_GENOME_INDEX_FILES) $(STAR_NOVEL_SPLICE_FILES)
index: hisat_index star_index

hisat_align: $(HISAT_BAM_FILES) $(HISAT_BAI_FILES)
star_align: $(STAR_BAM_FILES) $(STAR_BAI_FILES)
align: hisat_align star_align

hisat_qc: $(HISAT_QC_FILES)
hisat_multiqc: $(HISAT_MULTIQC_FILES)
star_qc: $(STAR_QC_FILES)
star_multiqc: $(STAR_MULTIQC_FILES)
qc: hisat_qc star_qc
multiqc: hisat_multiqc star_multiqc

hisat_count: $(HISAT_COUNT_FILES)
star_count: $(STAR_COUNT_FILES)
count: hisat_count star_count

hisat_coverage: $(HISAT_COVERAGE_FILES)
star_coverage: $(STAR_COVERAGE_FILES)
coverage: hisat_coverage star_coverage

hisat_cufflinks_ass: $(HISAT_CUFFLINKS_ASS_MERGED)
star_cufflinks_ass: $(STAR_CUFFLINKS_ASS_MERGED)
cufflinks_ass: star_cufflinks_ass hisat_cufflinks_ass

hisat_cufflinks_count: $(HISAT_CUFF_COUNT_FILES)
star_cufflinks_count: $(STAR_CUFF_COUNT_FILES)
cufflinks_count: hisat_cufflinks_count star_cufflinks_count

star_cufflinks_quant: $(ORIG_STAR_CUFF_QUANT_FILES) $(NEW_STAR_CUFF_QUANT_FILES)
cufflinks_quant: star_cufflinks_quant

star_cufflinks_norm: $(ORIG_STAR_CUFF_NORM_FILES)
cufflinks_norm: star_cufflinks_norm

# 01 - build index
$(HISAT_GENOME_INDEX_FILES): $(GENOME_FILE)
	@mkdir -p $(dir $@)
	$(HISAT_DOCKER) hisat2-build -f $< $(HISAT_GENOME_INDEX_DIR)/$(PREFIX)

$(HISAT_SPLICE_FILE): $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	$(HISAT_DOCKER) extract_splice_sites.py $< > $@

$(HISAT_NOVEL_SPLICE_FILE): $(HISAT_SPLICE_FILE) $(HISAT_GENOME_INDEX_FILES) $(READS1) $(READS2)
	@mkdir -p $(dir $@)
	$(HISAT_DOCKER) hisat2 \
	  --threads $(NCPU) \
	  --known-splicesite-infile $(HISAT_SPLICE_FILE) \
		--minins 0 \
		--maxins 220 \
		--min-intronlen 20 \
		--max-intronlen 10000 \
		--rna-strandness RF \
		--fr \
		--no-mixed \
		--no-discordant \
		--no-temp-splicesite \
		--dta-cufflinks \
		--novel-splicesite-outfile $@ \
		-x $(HISAT_GENOME_INDEX_DIR)/$(PREFIX) \
		-1 $(subst $(SPACE),$(COMMA),$(READS1)) \
		-2 $(subst $(SPACE),$(COMMA),$(READS2)) \
		-S /dev/null

$(STAR_GENOME_INDEX_DIR)/SA: $(GENOME_FILE) $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	$(STAR_DOCKER) STAR \
    --runThreadN $(NCPU) \
    --runMode genomeGenerate \
    --genomeDir $(STAR_GENOME_INDEX_DIR) \
    --genomeFastaFiles $(word 1, $^) \
    --sjdbGTFfile $(word 2, $^) \
    --sjdbOverhang 124 \
    --genomeSAindexNbases 11

$(STAR_GENOME_INDEX_DIR)/%.SJ.out.tab: $(DATA)/%$(READS1_PATTERN) $(DATA)/%$(READS2_PATTERN) $(STAR_GENOME_INDEX_FILES)
	@mkdir -p $(dir $@)
	$(STAR_DOCKER) STAR \
		--runThreadN $(NCPU) \
		--readFilesCommand zcat \
		--genomeDir $(STAR_GENOME_INDEX_DIR) \
		--outSAMtype None \
		--outSAMmode None \
		--outSJfilterReads All \
		--outSJfilterOverhangMin 5 1 1 1 \
		--outSJfilterCountUniqueMin 10 5 5 5 \
		--outSJfilterDistToOtherSJmin 5 0 0 3 \
		--outSJfilterIntronMaxVsReadN 0 1 500 5000 10000 20000 \
		--alignIntronMin 5 \
		--alignIntronMax 10000 \
		--alignSJoverhangMin 1 \
		--alignSJDBoverhangMin 1 \
		--alignSoftClipAtReferenceEnds No \
		--outFileNamePrefix $(basename $(basename $(basename $@))). \
		--readFilesIn $(word 1, $^) $(word 2, $^)

# 02 - align to reference
#		--no-temp-splicesite
$(HISAT_BAM_COMPLETE_TARGETS): $(DATA)/%$(READS1_PATTERN) $(DATA)/%$(READS2_PATTERN) $(HISAT_SPLICE_FILE) $(HISAT_NOVEL_SPLICE_FILE) $(HISAT_GENOME_INDEX_FILES)
	@mkdir -p $(dir $@)
	$(HISAT_DOCKER) hisat2 \
	  --threads $(NCPU) \
	  --known-splicesite-infile $(HISAT_SPLICE_FILE) \
		--novel-splicesite-infile $(HISAT_NOVEL_SPLICE_FILE) \
		--minins 0 \
		--maxins 220 \
		--min-intronlen 20 \
		--max-intronlen 10000 \
		--rna-strandness FR \
		--fr \
		--no-mixed \
		--no-discordant \
		-x $(HISAT_GENOME_INDEX_DIR)/$(PREFIX) \
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

$(HISAT_ALIGN_DIR)/%.bam: $(HISAT_ALIGN_DIR)/%.complete.bam
	$(SAMTOOLS_DOCKER) samtools view -bh -F 0x0004 -@ $(NCPU) $< > $@.tmp \
		&& mv $@.tmp $@

$(HISAT_ALIGN_DIR)/%.bai: $(HISAT_ALIGN_DIR)/%
	@mkdir -p $(dir $@)
	$(SAMTOOLS_DOCKER) samtools index $<


$(STAR_ALIGN_DIR)/%.bam: $(DATA)/%$(READS1_PATTERN) $(DATA)/%$(READS2_PATTERN) $(STAR_GENOME_INDEX_FILES) $(STAR_GENOME_INDEX_DIR)/*.SJ.out.tab
	@mkdir -p $(dir $@)
	$(STAR_DOCKER) STAR \
		--runThreadN $(NCPU) \
		--readFilesCommand zcat \
		--genomeDir $(STAR_GENOME_INDEX_DIR) \
		--sjdbFileChrStartEnd $(STAR_GENOME_INDEX_DIR)/*SJ.out.tab \
		--outSAMtype BAM Unsorted \
		--outBAMcompression 1 \
		--outSJfilterReads All \
		--outSJfilterCountUniqueMin 10 5 5 5 \
		--outSJfilterIntronMaxVsReadN 0 1 500 5000 10000 20000 \
		--alignIntronMin 5 \
		--alignIntronMax 10000 \
		--alignSJoverhangMin 10 \
		--alignSJDBoverhangMin 1 \
		--alignSoftClipAtReferenceEnds No \
		--outFilterType BySJout \
		--outFilterMultimapNmax 1 \
		--outFilterMismatchNmax 10 \
		--outFilterMismatchNoverLmax 0.2 \
		--outMultimapperOrder Random \
		--outSAMattributes All \
		--outSAMstrandField intronMotif\
		--outSAMattrIHstart 0 \
		--outSAMmapqUnique 50 \
		--outFileNamePrefix $(basename $@). \
		--readFilesIn $(word 1, $^) $(word 2, $^) \
	&& $(SAMTOOLS_DOCKER) samtools view \
	  	-uT $(GENOME_FILE) \
			$(basename $@).Aligned.out.bam \
		| $(SAMTOOLS_DOCKER) samtools sort \
		  -O BAM \
			-@ $(NCPU) \
			-l 9 \
			-o $@ \
	&& rm $(basename $@).Aligned.out.bam

$(STAR_ALIGN_DIR)/%.bai: $(STAR_ALIGN_DIR)/%
	@mkdir -p $(dir $@)
	$(SAMTOOLS_DOCKER) samtools index $<

# 03 - QC

$(HISAT_QC_DIR)/%.bam.bam_stat: $(HISAT_ALIGN_DIR)/%.complete.bam $(HISAT_ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) bam_stat.py -i $< > $@.tmp \
		&& mv $@.tmp $@

$(HISAT_CLIPPING_PROFILE_TARGETS): $(HISAT_ALIGN_DIR)/%.complete.bam $(HISAT_ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) clipping_profile.py \
		-i $< \
		-s "PE" \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(HISAT_DELETION_PROFILE_TARGETS): $(HISAT_ALIGN_DIR)/%.complete.bam $(HISAT_ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) deletion_profile.py \
		-i $< \
		-l 125 \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(HISAT_INNER_DISTANCE_TARGETS): $(HISAT_ALIGN_DIR)/%.complete.bam $(HISAT_ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) inner_distance.py \
		-i $< \
		--refgene=$(ANNOTATION_BED) \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(HISAT_INSERTION_PROFILE_TARGETS): $(HISAT_ALIGN_DIR)/%.complete.bam $(HISAT_ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) insertion_profile.py \
		-i $< \
		-s PE \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(HISAT_JUNCTION_ANNOTATION_TARGETS): $(HISAT_ALIGN_DIR)/%.complete.bam $(HISAT_ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) junction_annotation.py \
		-i $< --refgene=$(ANNOTATION_BED) \
		--min-intron=5 \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(HISAT_JUNCTION_SATURATION_TARGETS): $(HISAT_ALIGN_DIR)/%.complete.bam $(HISAT_ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) junction_saturation.py \
		-i $< \
		--refgene=$(ANNOTATION_BED) \
		--min-intron=5 \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(HISAT_MISMATCH_PROFILE_TARGETS): $(HISAT_ALIGN_DIR)/%.complete.bam $(HISAT_ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) mismatch_profile.py \
		-i $< \
		-l 125 \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(HISAT_READ_DUPLICATION_TARGETS): $(HISAT_ALIGN_DIR)/%.complete.bam $(HISAT_ALIGN_DIR)/%.complete.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) read_duplication.py \
		-i $< \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(HISAT_GENEBODY_COVERAGE_TARGETS): $(HISAT_BAM_COMPLETE_FILES) $(HISAT_BAI_COMPLETE_FILES)
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) geneBody_coverage.py \
		-i $(subst $(SPACE),$(COMMA),$(HISAT_BAM_COMPLETE_FILES)) \
		--refgene=$(ANNOTATION_BED) \
		-o $(HISAT_QC_DIR)/combined

$(HISAT_QC_DIR)/multiqc_report.html: $(HISAT_BAM_COMPLETE_FILES) $(HISAT_QC_FILES)
	@mkdir -p $(dir $@)
	$(MULTIQC_DOCKER) multiqc -o $(dir $@) $(HISAT_QC_DIR) $(HISAT_ALIGN_DIR)

$(STAR_QC_DIR)/%.bam.bam_stat: $(STAR_ALIGN_DIR)/%.bam $(STAR_ALIGN_DIR)/%.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) bam_stat.py -i $< > $@.tmp \
		&& mv $@.tmp $@

$(STAR_CLIPPING_PROFILE_TARGETS): $(STAR_ALIGN_DIR)/%.bam $(STAR_ALIGN_DIR)/%.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) clipping_profile.py \
		-i $< \
		-s "PE" \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(STAR_DELETION_PROFILE_TARGETS): $(STAR_ALIGN_DIR)/%.bam $(STAR_ALIGN_DIR)/%.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) deletion_profile.py \
		-i $< \

		-l 125 \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(STAR_INNER_DISTANCE_TARGETS): $(STAR_ALIGN_DIR)/%.bam $(STAR_ALIGN_DIR)/%.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) inner_distance.py \
		-i $< \
		--refgene=$(ANNOTATION_BED) \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(STAR_INSERTION_PROFILE_TARGETS): $(STAR_ALIGN_DIR)/%.bam $(STAR_ALIGN_DIR)/%.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) insertion_profile.py \
		-i $< \
		-s PE \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(STAR_JUNCTION_ANNOTATION_TARGETS): $(STAR_ALIGN_DIR)/%.bam $(STAR_ALIGN_DIR)/%.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) junction_annotation.py \
		-i $< --refgene=$(ANNOTATION_BED) \
		--min-intron=5 \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(STAR_JUNCTION_SATURATION_TARGETS): $(STAR_ALIGN_DIR)/%.bam $(STAR_ALIGN_DIR)/%.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) junction_saturation.py \
		-i $< \
		--refgene=$(ANNOTATION_BED) \
		--min-intron=5 \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(STAR_MISMATCH_PROFILE_TARGETS): $(STAR_ALIGN_DIR)/%.bam $(STAR_ALIGN_DIR)/%.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) mismatch_profile.py \
		-i $< \
		-l 125 \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(STAR_READ_DUPLICATION_TARGETS): $(STAR_ALIGN_DIR)/%.bam $(STAR_ALIGN_DIR)/%.bam.bai
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) read_duplication.py \
		-i $< \
		-o $(addprefix $(dir $@), $(notdir $(basename $(basename $@))))

$(STAR_GENEBODY_COVERAGE_TARGETS): $(STAR_BAM_FILES) $(STAR_BAI_FILES)
	@mkdir -p $(dir $@)
	$(RSEQC_DOCKER) geneBody_coverage.py \
		-i $(subst $(SPACE),$(COMMA),$(STAR_BAM_FILES)) \
		--refgene=$(ANNOTATION_BED) \
		-o $(STAR_QC_DIR)/combined

$(STAR_QC_DIR)/multiqc_report.html: $(STAR_BAM_FILES) $(STAR_QC_FILES)
	@mkdir -p $(dir $@)
	rm -rf $(STAR_QC_DIR)/multiqc*
	$(MULTIQC_DOCKER) multiqc -o $(dir $@) $(STAR_QC_DIR) $(STAR_ALIGN_DIR)

# 04 - count reads

# Note, because we sorted by position htseq-count uses a huge amount of RAM to count fragments.
# Inplanta counts used about 8 GB of RAM to count, invitro counts used 15 GB of RAM, if this is too much for your computer i suggest looking at a different program or temporarily sorting by name.
# Also by default, htseq-count sets a maximum ram useage of ~3GB, you need to edit the package following http://seqanswers.com/forums/showthread.php?p=197997 .

$(HISAT_COUNT_FILES): $(HISAT_BAM_FILES) $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	$(SUBREAD_DOCKER) featureCounts -a $(ANNOTATION_FILE) -t exon -s 2 -p -B -C -T $(NCPU) -o $@ $(HISAT_BAM_FILES)

$(STAR_COUNT_FILES): $(STAR_BAM_FILES) $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	$(SUBREAD_DOCKER) featureCounts -a $(ANNOTATION_FILE) -t exon -s 2 -p -B -C -T $(NCPU) -o $@ $(STAR_BAM_FILES)

# 05 - Coverage

$(HISAT_COVERAGE_DIR)/%-norm_factor.txt: $(HISAT_ALIGN_DIR)/%.bam
	@mkdir -p $(dir $@)
	echo "10^6/$$(samtools view -F 0x40 $< | cut -f1 | sort | uniq | wc -l)" | bc -l > $@.tmp && \
		mv $@.tmp $@

$(HISAT_COVERAGE_DIR)/%-coverage.bedgraph: $(HISAT_ALIGN_DIR)/%.bam $(GENOME_FILE) $(HISAT_COVERAGE_DIR)/%-norm_factor.txt
	@mkdir -p $(dir $@)
	bedtools genomecov -bga -split -trackline -scale $(shell cat $(word 3, $^)) -ibam $(word 1, $^) -g $(word 2, $^) > $@.tmp \
	  && mv $@.tmp $@

$(HISAT_COVERAGE_DIR)/%-coverage-fwd.bedgraph: $(HISAT_ALIGN_DIR)/%.bam $(GENOME_FILE) $(HISAT_COVERAGE_DIR)/%-norm_factor.txt
	@mkdir -p $(dir $@)
	bedtools genomecov -bga -split -trackline -scale $(shell cat $(word 3, $^)) -strand "+" -ibam $(word 1, $^) -g $(word 2, $^) > $@.tmp \
	  && mv $@.tmp $@

$(HISAT_COVERAGE_DIR)/%-coverage-rev.bedgraph: $(HISAT_ALIGN_DIR)/%.bam $(GENOME_FILE) $(HISAT_COVERAGE_DIR)/%-norm_factor.txt
	@mkdir -p $(dir $@)
	bedtools genomecov -bga -split -trackline -scale $(shell cat $(word 3, $^)) -strand "-" -ibam $(word 1, $^) -g $(word 2, $^) > $@.tmp \
	  && mv $@.tmp $@


$(STAR_COVERAGE_DIR)/%-norm_factor.txt: $(STAR_ALIGN_DIR)/%.bam
	@mkdir -p $(dir $@)
	echo "10^6/$$(samtools view -F 0x40 $< | cut -f1 | sort | uniq | wc -l)" | bc -l > $@.tmp && \
		mv $@.tmp $@

$(STAR_COVERAGE_DIR)/%-coverage.bedgraph: $(STAR_ALIGN_DIR)/%.bam $(GENOME_FILE) $(STAR_COVERAGE_DIR)/%-norm_factor.txt
	@mkdir -p $(dir $@)
	bedtools genomecov -bga -split -trackline -scale $(shell cat $(word 3, $^)) -ibam $(word 1, $^) -g $(word 2, $^) > $@.tmp \
	  && mv $@.tmp $@

$(STAR_COVERAGE_DIR)/%-coverage-fwd.bedgraph: $(STAR_ALIGN_DIR)/%.bam $(GENOME_FILE) $(STAR_COVERAGE_DIR)/%-norm_factor.txt
	@mkdir -p $(dir $@)
	bedtools genomecov -bga -split -trackline -scale $(shell cat $(word 3, $^)) -strand "+" -ibam $(word 1, $^) -g $(word 2, $^) > $@.tmp \
	  && mv $@.tmp $@

$(STAR_COVERAGE_DIR)/%-coverage-rev.bedgraph: $(STAR_ALIGN_DIR)/%.bam $(GENOME_FILE) $(STAR_COVERAGE_DIR)/%-norm_factor.txt
	@mkdir -p $(dir $@)
	bedtools genomecov -bga -split -trackline -scale $(shell cat $(word 3, $^)) -strand "-" -ibam $(word 1, $^) -g $(word 2, $^) > $@.tmp \
	  && mv $@.tmp $@

# 06a - cufflinks assemble

$(HISAT_CUFFLINKS_ASS_DIR)/%/transcripts.gtf: $(HISAT_ALIGN_DIR)/%.bam $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	cufflinks \
	  -p $(NCPU) \
		-g $(ANNOTATION_FILE) \
		-b $(GENOME_FILE) \
		--library-type=fr-firststrand \
		--max-intron-length 10000 \
		--min-intron-length 5 \
		-o $(dir $@) \
	  $(word 1, $^)

$(STAR_CUFFLINKS_ASS_DIR)/%/transcripts.gtf: $(STAR_ALIGN_DIR)/%.bam $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	cufflinks \
	  -p $(NCPU) \
		-g $(ANNOTATION_FILE) \
		-b $(GENOME_FILE) \
		--library-type=fr-firststrand \
		--max-intron-length 10000 \
		--min-intron-length 5 \
		-o $(dir $@) \
	  $(word 1, $^)


# 06b - cufflinks merge

$(HISAT_CUFFLINKS_ASS_MERGED): $(HISAT_CUFFLINKS_ASS_FILES) $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)

	@echo "$(HISAT_CUFFLINKS_ASS_FILES)" | tr " " "\n" > $(dir $@)/to_be_merged.txt
	cuffmerge \
		--merge
	  -p $(NCPU) \
		-g $(ANNOTATION_FILE) \
		-s $(GENOME_FILE) \
		-o $(dir $@) \
		$(dir $@)/to_be_merged.txt

$(STAR_CUFFLINKS_ASS_MERGED): $(STAR_CUFFLINKS_ASS_FILES) $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)

	@echo "$(STAR_CUFFLINKS_ASS_FILES)" | tr " " "\n" > $(dir $@)/to_be_merged.txt
	cuffmerge \
	  -p $(NCPU) \
    -g $(ANNOTATION_FILE) \
    -s $(GENOME_FILE) \
    -o $(dir $@) \
    $(dir $@)/to_be_merged.txt

# 07 - count cufflinks

$(HISAT_CUFF_COUNT_FILES): $(HISAT_BAM_FILES) $(HISAT_CUFFLINKS_ASS_MERGED)
	@mkdir -p $(dir $@)
	$(SUBREAD_DOCKER) featureCounts -a $(HISAT_CUFFLINKS_ASS_MERGED) -t exon -s 2 -p -B -C -T $(NCPU) -o $@ $(HISAT_BAM_FILES)

$(STAR_CUFF_COUNT_FILES): $(STAR_BAM_FILES) $(STAR_CUFFLINKS_ASS_MERGED)
	@mkdir -p $(dir $@)
	$(SUBREAD_DOCKER) featureCounts -a $(STAR_CUFFLINKS_ASS_MERGED) -t exon -s 2 -p -B -C -T $(NCPU) -o $@ $(STAR_BAM_FILES)


# 08

$(ORIG_STAR_CUFF_QUANT_DIR)/%/abundances.cxb: $(STAR_ALIGN_DIR)/%.bam $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	cuffquant \
	  -p $(NCPU) \
	  --frag-bias-correct $(GENOME_FILE) \
		--multi-read-correct \
		--library-type=fr-firststrand \
		--output-dir $(dir $@) \
		$(ANNOTATION_FILE) \
	  $(word 1, $^)

$(NEW_STAR_CUFF_QUANT_DIR)/%/abundances.cxb: $(STAR_ALIGN_DIR)/%.bam $(STAR_CUFFLINKS_ASS_MERGED)
	@mkdir -p $(dir $@)
	cuffquant \
	  -p $(NCPU) \
		--frag-bias-correct $(GENOME_FILE) \
		--multi-read-correct \
		--library-type=fr-firststrand \
		--output-dir $(dir $@) \
		$(STAR_CUFFLINKS_ASS_MERGED) \
	  $(word 1, $^)

$(ORIG_STAR_CUFF_NORM_FILES): $(ANNOTATION_FILE) $(ORIG_STAR_CUFF_QUANT_FILES)
	@mkdir -p $(dir $@)
	cuffnorm \
		-p $(NCPU) \
		--library-type=fr-firststrand \
		--library-norm-method=geometric \
		--output-dir $(dir $@) \
		$(ANNOTATION_FILE) \
		cufflinks_quant_orig/star/inplanta-3KO-020/abundances.cxb,cufflinks_quant_orig/star/inplanta-3KO-021/abundances.cxb,cufflinks_quant_orig/star/inplanta-3KO-022/abundances.cxb,cufflinks_quant_orig/star/inplanta-3KO-023/abundances.cxb \
		cufflinks_quant_orig/star/inplanta-Pf2-024/abundances.cxb,cufflinks_quant_orig/star/inplanta-Pf2-025/abundances.cxb,cufflinks_quant_orig/star/inplanta-Pf2-026/abundances.cxb,cufflinks_quant_orig/star/inplanta-Pf2-027/abundances.cxb \
		cufflinks_quant_orig/star/inplanta-SN15-016/abundances.cxb,cufflinks_quant_orig/star/inplanta-SN15-017/abundances.cxb,cufflinks_quant_orig/star/inplanta-SN15-018/abundances.cxb,cufflinks_quant_orig/star/inplanta-SN15-019/abundances.cxb \
		cufflinks_quant_orig/star/invitro-3KO-004/abundances.cxb,cufflinks_quant_orig/star/invitro-3KO-005/abundances.cxb,cufflinks_quant_orig/star/invitro-3KO-006/abundances.cxb,cufflinks_quant_orig/star/invitro-3KO-007/abundances.cxb \
		cufflinks_quant_orig/star/invitro-Pf2-008/abundances.cxb,cufflinks_quant_orig/star/invitro-Pf2-009/abundances.cxb,cufflinks_quant_orig/star/invitro-Pf2-010/abundances.cxb,cufflinks_quant_orig/star/invitro-Pf2-011/abundances.cxb \
		cufflinks_quant_orig/star/invitro-SN15-000/abundances.cxb,cufflinks_quant_orig/star/invitro-SN15-001/abundances.cxb,cufflinks_quant_orig/star/invitro-SN15-002/abundances.cxb,cufflinks_quant_orig/star/invitro-SN15-003/abundances.cxb



$(NEW_STAR_CUFF_NORM_FILES): $(ANNOTATION_FILE) $(NEW_STAR_CUFF_QUANT_FILES)
	@mkdir -p $(dir $@)
	cuffnorm \
		-p $(NCPU) \
		--library-type=fr-firststrand \
		--library-norm-method=geometric \
		--output-dir $(dir $@) \
		$(ANNOTATION_FILE) \
		cufflinks_quant_new/star/inplanta-3KO-020/abundances.cxb,cufflinks_quant_new/star/inplanta-3KO-021/abundances.cxb,cufflinks_quant_new/star/inplanta-3KO-022/abundances.cxb,cufflinks_quant_new/star/inplanta-3KO-023/abundances.cxb \
		cufflinks_quant_new/star/inplanta-Pf2-024/abundances.cxb,cufflinks_quant_new/star/inplanta-Pf2-025/abundances.cxb,cufflinks_quant_new/star/inplanta-Pf2-026/abundances.cxb,cufflinks_quant_new/star/inplanta-Pf2-027/abundances.cxb \
		cufflinks_quant_new/star/inplanta-SN15-016/abundances.cxb,cufflinks_quant_new/star/inplanta-SN15-017/abundances.cxb,cufflinks_quant_new/star/inplanta-SN15-018/abundances.cxb,cufflinks_quant_new/star/inplanta-SN15-019/abundances.cxb \
		cufflinks_quant_new/star/invitro-3KO-004/abundances.cxb,cufflinks_quant_new/star/invitro-3KO-005/abundances.cxb,cufflinks_quant_new/star/invitro-3KO-006/abundances.cxb,cufflinks_quant_new/star/invitro-3KO-007/abundances.cxb \
		cufflinks_quant_new/star/invitro-Pf2-008/abundances.cxb,cufflinks_quant_new/star/invitro-Pf2-009/abundances.cxb,cufflinks_quant_new/star/invitro-Pf2-010/abundances.cxb,cufflinks_quant_new/star/invitro-Pf2-011/abundances.cxb \
		cufflinks_quant_new/star/invitro-SN15-000/abundances.cxb,cufflinks_quant_new/star/invitro-SN15-001/abundances.cxb,cufflinks_quant_new/star/invitro-SN15-002/abundances.cxb,cufflinks_quant_new/star/invitro-SN15-003/abundances.cxb
