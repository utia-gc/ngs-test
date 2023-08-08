#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

#######################################
# Aliases
#######################################
function samtools() {
    singularity exec https://depot.galaxyproject.org/singularity/samtools%3A1.17--hd87286a_1 samtools "${@}"
}
export -f samtools

#######################################
# Constants
#######################################
readonly SUBSAMPLE_SEED="${1}"
readonly SUBSAMPLE_FRACTION="${2}"
readonly BAM="${3}"
readonly REGION="${4}"

#######################################
# Functions
#######################################

#######################################
# Downsample BAM file.
# Arguments:
#   Subsample seed, an integer, the seed to set for subsampling.
#   Subsample fraction, a float between 0 and 1, the fraction of reads/pairs to keep.
#   BAM file, a path.
#   Region, a string RNAME[:STARTPOS[-ENDPOS]], the region within which to consider alignments.
# Outputs:
#   Writes single read name each line.
# Returns:
#   0 if completed, non-zero on error.
#######################################
function downsample_bam() {
    samtools view \
        --subsample-seed "${1}" \
        --subsample "${2}" \
        "${3}" \
        "${4}"
}

#######################################
# Extract read names from a SAM file.
# Arguments:
#   SAM file.
# Outputs:
#   Writes single read name each line.
# Returns:
#   0 if completed, non-zero on error.
#######################################
function extract_sam_read_names() {
    cut -f 1 "${1}"
}

#######################################
# Main file.
# Globals:
#   Subsample seed, an integer, the seed to set for subsampling.
#   Subsample fraction, a float between 0 and 1, the fraction of reads/pairs to keep.
#   BAM file, a path.
#   Region, a string RNAME[:STARTPOS[-ENDPOS]], the region within which to consider alignments.
# Outputs:
#   Writes single read name each line.
# Returns:
#   0 if completed, non-zero on error.
#######################################
function main() {
    downsample_bam "${SUBSAMPLE_SEED}" "${SUBSAMPLE_FRACTION}" "${BAM}" "${REGION}" \
    | extract_sam_read_names -
}

main
