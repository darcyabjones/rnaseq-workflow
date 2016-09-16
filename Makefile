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
STAR_DOCKER=sudo docker run --rm -v $PWD:/data:z darcyabjones/star
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


# 02 - align to reference
# Complete ones have all reads (even unaligned, use this for QC), regular ones just have the mapped ones.

ALIGN_DIR=align
HISAT_ALIGN_DIR=$(ALIGN_DIR)/hisat2
HISAT_BAM_COMPLETE_TARGETS=$(addprefix $(HISAT_ALIGN_DIR)/, $(addsuffix .complete.bam, %)))
HISAT_BAM_COMPLETE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_BAM_COMPLETE_TARGETS)))
HISAT_BAI_COMPLETE_FILES=$(addsuffix .bai, $(HISAT_BAM_COMPLETE_FILES))

HISAT_BAM_FILES=$(addprefix $(HISAT_ALIGN_DIR)/, $(addsuffix .bam, $(SAMPLE_NAMES)))
HISAT_BAI_FILES=$(addsuffix .bai, $(HISAT_BAM_FILES))

# 03 - QC

QC_DIR=bam_qc

HISAT_QC_DIR=$(QC_DIR)/hisat2
HISAT_BAM_STAT_FILES=$(addprefix $(HISAT_QC_DIR)/, $(addsuffix .bam_stat, $(notdir $(HISAT_BAM_FILES))))

CLIPPING_PROFILE_EXTS=-clipping_profile.r .clipping_profile.R1.pdf .clipping_profile.R2.pdf .clipping_profile.xls
HISAT_CLIPPING_PROFILE_TARGETS=$(foreach e, $(CLIPPING_PROFILE_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_CLIPPING_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_CLIPPING_PROFILE_TARGETS)))

DELETION_PROFILE_EXTS=.deletion_profile.r .deletion_profile.txt
HISAT_DELETION_PROFILE_TARGETS=$(foreach e, $(DELETION_PROFILE_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_DELETION_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_DELETION_PROFILE_TARGETS)))

INNER_DISTANCE_EXTS=.inner_distance_freq.txt .inner_distance_plot.pdf .inner_distance_plot.r .inner_distance_plot.r
HISAT_INNER_DISTANCE_TARGETS=$(foreach e, $(INNER_DISTANCE_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_INNER_DISTANCE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_INNER_DISTANCE_TARGETS)))

INSERTION_PROFILE_EXTS=.insertion_profile.r .insertion_profile.R1.pdf .insertion_profile.R2.pdf .insertion_profile.xls
HISAT_INSERTION_PROFILE_TARGETS=$(foreach e, $(INSERTION_PROFILE_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_INSERTION_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_INSERTION_PROFILE_TARGETS)))

JUNCTION_ANNOTATION_EXTS=.junction.bed .junction_plot.r .junction.xls
HISAT_JUNCTION_ANNOTATION_TARGETS=$(foreach e, $(JUNCTION_ANNOTATION_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_JUNCTION_ANNOTATION_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_JUNCTION_ANNOTATION_TARGETS)))

JUNCTION_SATURATION_EXTS=.junctionSaturation_plot.pdf .junctionSaturation_plot.r
HISAT_JUNCTION_SATURATION_TARGETS=$(foreach e, $(JUNCTION_SATURATION_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_JUNCTION_SATURATION_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_JUNCTION_SATURATION_TARGETS)))

MISMATCH_PROFILE_EXTS=.mismatch_profile.pdf .mismatch_profile.r .mismatch_profile.xls
HISAT_MISMATCH_PROFILE_TARGETS=$(foreach e, $(MISMATCH_PROFILE_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_MISMATCH_PROFILE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_MISMATCH_PROFILE_TARGETS)))

READ_DUPLICATION_EXTS=.DupRate_plot.pdf .DupRate_plot.r .pos.DupRate.xls .seq.DupRate.xls
HISAT_READ_DUPLICATION_TARGETS=$(foreach e, $(READ_DUPLICATION_EXTS), $(addprefix $(HISAT_QC_DIR)/, $(addsuffix $(e), %)))
HISAT_READ_DUPLICATION_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_READ_DUPLICATION_TARGETS)))

GENEBODY_COVERAGE_EXTS=.geneBodyCoverage.curves.pdf .geneBodyCoverage.r .geneBodyCoverage.txt
HISAT_GENEBODY_COVERAGE_TARGETS=$(foreach e, $(GENEBODY_COVERAGE_EXTS), $(addsuffix $(e), $(HISAT_QC_DIR)/%))
HISAT_GENEBODY_COVERAGE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_GENEBODY_COVERAGE_TARGETS)))

HISAT_MULTIQC_FILES=$(HISAT_QC_DIR)/multiqc_report.html

HISAT_QC_FILES=$(HISAT_BAM_STAT_FILES) $(HISAT_CLIPPING_PROFILE_FILES) \
  $(HISAT_DELETION_PROFILE_FILES) $(HISAT_INNER_DISTANCE_FILES) \
	$(HISAT_INSERTION_PROFILE_FILES) $(HISAT_JUNCTION_ANNOTATION_FILES) \
	$(HISAT_JUNCTION_SATURATION_FILES) $(HISAT_MISMATCH_PROFILE_FILES) \
	$(HISAT_READ_DUPLICATION_FILES) $(HISAT_GENEBODY_COVERAGE_FILE)



# 04 - Count reads

COUNT_DIR=counts

HISAT_COUNT_DIR=$(COUNT_DIR)/hisat2
HISAT_COUNT_FILES=$(addprefix $(HISAT_COUNT_DIR)/, $(addsuffix -counts.tsv, $(SAMPLE_NAMES)))

# 05 - Calculate coverage

COVERAGE_DIR=coverage

HISAT_COVERAGE_DIR=$(COVERAGE_DIR)/hisat2
HISAT_COVERAGE_FILES=$(addprefix $(HISAT_COVERAGE_DIR)/, $(addsuffix -coverage.bedgraph, $(SAMPLE_NAMES)))

# 06 - Assemble transcripts.

STRINGTIE_DIR=stringtie
STRINGTIE_EXTS=.stringtie.gtf .stringtie.ctab

HISAT_STRINGTIE_DIR=$(STRINGTIE_DIR)/hisat2
HISAT_STRINGTIE_TARGETS=$(foreach e, $(ASSEMBLE_EXTS), $(addprefix $(HISAT_STRINGTIE_DIR)/, $(addsuffix $(e), %)))
HISAT_STRINGTIE_FILES=$(foreach s, $(SAMPLE_NAMES), $(subst %,$(s), $(HISAT_STRINGTIE_TARGETS)))


# Do the actual work
phony:
	$(HISAT_DOCKER) ls -l data

cleanall: all clean
	@rm -rf $^

clean: $(HISAT_BAM_COMPLETE_FILES) $(HISAT_BAI_COMPLETE_FILES)
	@rm -rf $^

all: index align qc count coverage

hisat_index: $(HISAT_GENOME_INDEX_FILES) $(HISAT_SPLICE_FILE) $(HISAT_NOVEL_SPLICE_FILE)
star_index:
index: hisat_index star_index

hisat_align: $(HISAT_BAM_FILES) $(HISAT_BAI_FILES)
star_align:
align: hisat_align star_align

hisat_qc: $(HISAT_QC_FILES) $(HISAT_MULTIQC_FILES)
star_qc:
qc: hisat_qc star_qc

hisat_count: $(HISAT_COUNT_FILES)
star_count:
count: hisat_count star_count

hisat_coverage: $(HISAT_COVERAGE_FILES)
star_coverage:
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
		--downstream-transcriptome-assembly \
		--novel-splicesite-outfile $@ \
		-x $(HISAT_GENOME_INDEX_DIR)/$(PREFIX) \
		-1 $(subst $(SPACE),$(COMMA),$(READS1)) \
		-2 $(subst $(SPACE),$(COMMA),$(READS2)) \
		-S /dev/null

#STAR \
#--runThreadN $(grep -c '^processor' /proc/cpuinfo) \
#--runMode genomeGenerate \
#--genomeDir genome_inde/star \
#--genomeFastaFiles data/Parastagonospora_nodorum_SN15_scaffolds.fasta \
#--sjdbGTFfile data/Parastagonospora_nodorum_SN15.gtf \
#--sjdbOverhang 124 \
#--genomeSAindexNbases 11


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
		--downstream-transcriptome-assembly \
		--rna-strandness FR \
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

#STAR --runThreadN $(grep -c '^processor' /proc/cpuinfo) --alignIntronMin 5 --genomeDir 03-index_genome/star --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $this_dir/Pn-ip-0${i}. --outSAMattributes All --readFilesIn 02-simulate_reads/inplanta_reads/sample_0${i}_1.fasta 02-simulate_reads/inplanta_reads/sample_0${i}_2.fasta
#STAR --runThreadN $(grep -c '^processor' /proc/cpuinfo) --sjdbFileChrStartEnd 04-align_reads/star-default/*SJ.out.tab --genomeDir 03-index_genome/star --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $this_dir/Pn-ip-0${i}. --outSAMattributes All --readFilesIn 02-simulate_reads/inplanta_reads/sample_0${i}_1.fasta 02-simulate_reads/inplanta_reads/sample_0${i}_2.fasta

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
		-i $(subst $(SPACE),$(COMMA),$(HISAT_BAM_COMPLETE_FILES))
		--refgene=$(ANNOTATION_BED) \
		-o $(HISAT_QC_DIR)/combined

$(HISAT_QC_DIR)/multiqc_report.html: $(HISAT_BAM_COMPLETE_FILES) $(HISAT_QC_FILES)
	@mkdir -p $(dir $@)
	$(MULTIQC_DOCKER) multiqc -o $(dir $@) $(HISAT_QC_DIR)

# 04 - count reads

# Note, because we sorted by position htseq-count uses a huge amount of RAM to count fragments.
# Inplanta counts used about 8 GB of RAM to count, invitro counts used 15 GB of RAM, if this is too much for your computer i suggest looking at a different program or temporarily sorting by name.
# Also by default, htseq-count sets a maximum ram useage of ~3GB, you need to edit the package following http://seqanswers.com/forums/showthread.php?p=197997 .

$(HISAT_COUNT_DIR)/%-counts.tsv: $(HISAT_ALIGN_DIR)/%.bam $(ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	$(HTSEQ_DOCKER) htseq-count --format=bam --order=pos --type=exon $(word 1, $^) $(word 2, $^) > $@.tmp \
	  && mv $@.tmp $@

# 05 - Coverage

$(HISAT_COVERAGE_DIR)/%-coverage.bedgraph: $(HISAT_ALIGN_DIR)/%.bam $(GENOME_FILE)
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
