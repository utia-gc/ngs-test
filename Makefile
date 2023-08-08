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

# Build the most "up-to-date" rule
build: data/temp/alignments/.alignments.sentinel data/temp/reads/SC3pv3_GEX_Human_PBMC/SC3pv3_GEX_Human_PBMC_fastqs.tar
.PHONY: build

# Clean up temporary and sentinel files.
# WARNING: Running make clean will require everything to be rebuilt.
clean:
> rm -rf data/temp/alignments
> rm -rf data/temp/reads
.PHONY: clean

# Fetch alignments files
data/temp/alignments/.alignments.sentinel: data/urls/alignments_urls.txt
> mkdir -p $(@D)
> wget --no-clobber --directory-prefix=$(@D) --input-file=data/urls/alignments_urls.txt
> touch $@

# Fetch reads files
data/temp/reads/SC3pv3_GEX_Human_PBMC/SC3pv3_GEX_Human_PBMC_fastqs.tar: data/urls/reads_urls.txt
> mkdir -p $(@D)
> wget --no-clobber --directory-prefix=$(@D) --input-file=data/urls/reads_urls.txt
