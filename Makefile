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

# Build the most "up-to-date" rule
build: data/temp/reads/subset/$(region)/.reads_$(region).sentinel data/temp/references/.references.sentinel
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
