SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

ifeq ($(origin .RECIPEPREFIX), undefined)
  $(error This Make does not support .RECIPEPREFIX. Please use GNU Make 4.0 or later)
endif
.RECIPEPREFIX = >

region := "chr21"
refseq_chr := "NC_000021.9"

# Build the most "up-to-date" rule
build: data/reads/raw/.L001.sentinel data/reads/raw/.L002.sentinel data/references/GRCh38.p14_genomic_$(refseq_chr).gtf.gz data/references/GRCh38.p14_genomic_$(refseq_chr).fna.gz
.PHONY: build

# Clean up temporary and sentinel files.
# WARNING: Running make clean will require everything to be rebuilt.
clean:
> rm -rf data/temp/alignments
> rm -rf data/temp/reads
> rm -rf data/temp/references
.PHONY: clean

# Fetch alignments files
data/temp/alignments/.alignments.sentinel: data/urls/alignments_urls.txt
> mkdir -p $(@D)
> wget --no-clobber --directory-prefix=$(@D) --input-file=data/urls/alignments_urls.txt
> touch $@

# Fetch reads files
data/temp/reads/SC3pv3_GEX_Human_PBMC_fastqs.tar: data/urls/reads_urls.txt
> mkdir -p $(@D)
> wget --no-clobber --directory-prefix=$(@D) --input-file=data/urls/reads_urls.txt

# Fetch references files
data/temp/references/.references.sentinel: data/urls/references_urls.txt
> mkdir -p $(@D)
> wget --no-clobber --directory-prefix=$(@D) --input-file=$<
> touch $@

# Extract reads from archive
data/temp/reads/full/.extract_reads.sentinel: data/temp/reads/SC3pv3_GEX_Human_PBMC_fastqs.tar
> mkdir -p $(@D)
> tar -xvf $< --directory=$(@D)
> touch $@

# Get names of reads aligned to region
data/temp/alignments/$(region)_read_names.txt: data/temp/alignments/.alignments.sentinel
> mkdir -p $(@D)
> bash src/bash/get_bam_read_names.bash 0 1 data/temp/alignments/SC3pv3_GEX_Human_PBMC_possorted_genome_bam.bam $(region) > $@

# Subset fastq files by read names
fastqs := $(shell find data/temp/reads/full -name "*.fastq.gz")
data/temp/reads/subset/$(region)/.reads_$(region).sentinel: $(fastqs) data/temp/alignments/$(region)_read_names.txt
> mkdir -p $(@D)
> for fastq in $(fastqs); do
> 	outname=$$(bash src/bash/append_to_sample_name.bash "$${fastq}" $(region))
> 	bash src/bash/subset_fastq_by_read_names.bash "$${fastq}" data/temp/alignments/$(region)_read_names.txt > $(@D)/"$${outname}"
> done
> touch $@

# Sample spoofed lane 001 reads
subset_fastqs := $(shell find data/temp/reads/subset -name "*.fastq.gz")
data/reads/raw/.L001.sentinel: $(subset_fastqs)
> mkdir -p $(@D)
> nreads=100000
> for fastq in $(subset_fastqs); do
> 	outname=$$(bash src/bash/append_to_sample_name.bash "$${fastq}" "$${nreads}")
> 	outname=$$(bash src/bash/change_lane_number.bash "$${outname}" 001)
> 	bash src/bash/downsample_reads.bash 100 "$${fastq}" "$${nreads}" > $(@D)/"$${outname}"
> done
> touch $@

# Sample spoofed lane 002 reads
subset_fastqs := $(shell find data/temp/reads/subset -name "*.fastq.gz")
data/reads/raw/.L002.sentinel: $(subset_fastqs)
> mkdir -p $(@D)
> nreads=100000
> for fastq in $(subset_fastqs); do
> 	outname=$$(bash src/bash/append_to_sample_name.bash "$${fastq}" "$${nreads}")
> 	outname=$$(bash src/bash/change_lane_number.bash "$${outname}" 002)
> 	bash src/bash/downsample_reads.bash 2023 "$${fastq}" "$${nreads}" > $(@D)/"$${outname}"
> done
> touch $@

# Filter GTF for only features in chromosome of interest
data/references/GRCh38.p14_genomic_$(refseq_chr).gtf.gz: data/temp/references/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
> mkdir -p $(@D)
> outname=$(@D)/$$(basename $@ ".gz")
> zgrep "^\#gtf-version" $< > "$${outname}"
> zgrep "^\#!" $< >> "$${outname}"
> awk -F '\t' '$$1 == $(refseq_chr)' <(zcat $<) >> "$${outname}"
> echo "###" >> "$${outname}"
> gzip --force "$${outname}"
> gunzip --keep $@

# Filter FASTA for only features in chromosome of interest
data/references/GRCh38.p14_genomic_$(refseq_chr).fna.gz: data/temp/references/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
> mkdir -p $(@D)
> echo $(refseq_chr) > $(@D)/name.lst
> singularity exec https://depot.galaxyproject.org/singularity/seqtk%3A1.4--he4a0461_1 seqtk subseq $< $(@D)/name.lst | gzip --force --stdout > $@
> gunzip --keep $@
> rm -f $(@D)/name.lst
