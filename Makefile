C=$(shell pwd)
SPACE :=
SPACE +=
COMMA := ,
NCPU=4

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
#STAR_DOCKER=sudo docker run --rm -v $PWD:/data:z darcyabjones/star
STAR_DOCKER=
#SAMTOOLS_DOCKER=sudo docker run --rm -v $$PWD:/data:z darcyabjones/samtools
SAMTOOLS_DOCKER=
RSEQC_DOCKER=
MULTIQC_DOCKER=
HTSEQ_DOCKER=
#STRINGTIE_DOCKER=sudo docker run --rm -v $$PWD:/data:z darcyabjones/stringtie
STRINGTIE_DOCKER=

# 01 - build index

GENOME_INDEX_DIR=genome_index

HISAT_GENOME_INDEX_DIR=$(GENOME_INDEX_DIR)/hisat2
HISAT_GENOME_INDEX_FILES=$(foreach i, 1, $(HISAT_GENOME_INDEX_DIR)/$(PREFIX).$i.ht2)
HISAT_SPLICE_FILE=$(HISAT_GENOME_INDEX_DIR)/splice_sites.txt
HISAT_EXTRACT_SPLICE=$(HISAT_DOCKER) extract_splice_sites.py $(1) > $(HISAT_SPLICE_FILE)
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
HISAT_BAM_COMPLETE_TARGETS=$(addprefix $(HISAT_ALIGN_DIR)/, $(addsuffix .complete.bam, %)))
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
HISAT_COUNT_FILES=$(addprefix $(HISAT_COUNT_DIR)/, $(addsuffix -counts.tsv, $(SAMPLE_NAMES)))

STAR_COUNT_DIR=$(COUNT_DIR)/star
STAR_COUNT_FILES=$(addprefix $(STAR_COUNT_DIR)/, $(addsuffix -counts.tsv, $(SAMPLE_NAMES)))

# 05 - Calculate coverage

COVERAGE_DIR=coverage

HISAT_COVERAGE_DIR=$(COVERAGE_DIR)/hisat2
HISAT_COVERAGE_FILES=$(addprefix $(HISAT_COVERAGE_DIR)/, $(addsuffix -coverage.bedgraph, $(SAMPLE_NAMES)))

STAR_COVERAGE_DIR=$(COVERAGE_DIR)/star
STAR_COVERAGE_FILES=$(addprefix $(STAR_COVERAGE_DIR)/, $(addsuffix -coverage.bedgraph, $(SAMPLE_NAMES)))

# 06 - Assemble transcripts.

STRINGTIE_DIR=stringtie
STRINGTIE_EXTS=.stringtie.gtf .stringtie.ctab

HISAT_STRINGTIE_DIR=$(STRINGTIE_DIR)/hisat2
HISAT_STRINGTIE_TARGETS=$(foreach e, $(STRINGTIE_EXTS), $(addprefix $(HISAT_STRINGTIE_DIR)/, $(addsuffix $(e), %)))
HISAT_STRINGTIE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_STRINGTIE_TARGETS)))

STAR_STRINGTIE_DIR=$(STRINGTIE_DIR)/star
STAR_STRINGTIE_TARGETS=$(foreach e, $(STRINGTIE_EXTS), $(addprefix $(STAR_STRINGTIE_DIR)/, $(addsuffix $(e), %)))
STAR_STRINGTIE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(STAR_STRINGTIE_TARGETS)))

# Do the actual work
phony:
	@echo $(STAR_NOVEL_SPLICE_FILES)

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
multiqc: hisat_multiqc

hisat_count: $(HISAT_COUNT_FILES)
star_count: $(STAR_COUNT_FILES)
count: hisat_count star_count

hisat_coverage: $(HISAT_COVERAGE_FILES)
star_coverage: $(STAR_COVERAGE_FILES)
coverage: hisat_coverage star_coverage

hisat_stringtie: $(HISAT_STRINGTIE_FILES)
star_stringtie:
stringtie: hisat_stringtie star_stringtie

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
		--rna-strandness FR \
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
		--no-mixed \
		--no-discordant \
		--no-temp-splicesite \
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


$(STAR_ALIGN_DIR)%.bam: $(DATA)/%$(READS1_PATTERN) $(DATA)/%$(READS2_PATTERN) $(STAR_GENOME_INDEX_FILES) $(STAR_GENOME_INDEX_DIR)/*.SJ.out.tab
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

$(HISAT_COUNT_DIR)/%-counts.tsv: $(HISAT_ALIGN_DIR)/%.bam $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	$(HTSEQ_DOCKER) htseq-count --format=bam --order=pos --type=exon --stranded=yes --mode=union $(word 1, $^) $(word 2, $^) > $@.tmp \
	  && mv $@.tmp $@

$(STAR_COUNT_DIR)/%-counts.tsv: $(STAR_ALIGN_DIR)/%.bam $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	$(HTSEQ_DOCKER) htseq-count --format=bam --order=pos --type=exon --stranded=yes --mode=union $(word 1, $^) $(word 2, $^) > $@.tmp \
	  && mv $@.tmp $@


# 05 - Coverage

$(HISAT_COVERAGE_DIR)/%-coverage.bedgraph: $(HISAT_ALIGN_DIR)/%.bam $(GENOME_FILE)
	@mkdir -p $(dir $@)
	bedtools genomecov -bga -split -trackline -ibam $(word 1, $^) -g $(word 2, $^) > $@.tmp \
	  && mv $@.tmp $@

$(STAR_COVERAGE_DIR)/%-coverage.bedgraph: $(STAR_ALIGN_DIR)/%.bam $(GENOME_FILE)
	@mkdir -p $(dir $@)
	bedtools genomecov -bga -split -trackline -ibam $(word 1, $^) -g $(word 2, $^) > $@.tmp \
	  && mv $@.tmp $@

# 06 - stringtie

$(HISAT_STRINGTIE_TARGETS): $(HISAT_ALIGN_DIR)/%.bam $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	$(STRINGTIE_DOCKER) stringtie \
	  $(word 1, $^) \
	  -p $(NCPU) \
		-G $(ANNOTATION_FILE) \
		-e \
		-b $(dir $@)/$(notdir $(basename $(word 1, $^))) \
		-o $(dir $@)/$(notdir $(basename $(word 1, $^))).stringtie.gtf
