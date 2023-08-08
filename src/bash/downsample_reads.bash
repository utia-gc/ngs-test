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
readonly SEED="${1}"
readonly FASTQ="${2}"
readonly NREADS="${3}"

#######################################
# Functions
#######################################

#######################################
# Downsample fastq file
# Arguments:
#   Random seed RNG.
#   Input fastq path.
#   Number reads to sample.
# Outputs:
#   Writes gzipped fastq file to stdout.
# Returns:
#   0 if completed, non-zero on error.
#######################################
function downsample_fastq() {
    seqtk sample -s"${1}" "${2}" "${3}" \
    | gzip
}

#######################################
# Main program.
# Globals:
#   SEED
#   FASTQ
#   NREADS
# Outputs:
#   Writes gzipped fastq file to stdout.
# Returns:
#   0 if completed, non-zero on error.
#######################################
function main() {
    downsample_fastq "${SEED}" "${FASTQ}" "${NREADS}"
}

main
