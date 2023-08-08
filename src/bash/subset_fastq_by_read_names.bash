#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

#######################################
# Aliases
#######################################
function seqtk() {
    singularity exec https://depot.galaxyproject.org/singularity/seqtk%3A1.4--he4a0461_1 seqtk "${@}"
}
export -f seqtk

#######################################
# Constants
#######################################
readonly FASTQ="${1}"
readonly READ_NAMES="${2}"

#######################################
# Functions
#######################################

#######################################
# Subset fastq from read names.
# Arguments:
#   Input fastq path.
#   Input reads names file path.
# Outputs:
#   Writes gzipped fastq file to output fastq path.
# Returns:
#   0 if completed, non-zero on error.
#######################################
function subset_fastq_by_read_names() {
    seqtk subseq "${1}" "${2}" \
    | gzip
}

#######################################
# Globals:
#   Input fastq path.
#   Input reads names file path.
# Main program.
# Returns:
#   0 if completed, non-zero on error.
#######################################
function main() {
    subset_fastq_by_read_names "${FASTQ}" "${READ_NAMES}"
}

main
