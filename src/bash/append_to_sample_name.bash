#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

#######################################
# Constants
#######################################
readonly FASTQ="${1}"
readonly REGION="${2}"

#######################################
# Functions
#######################################

#######################################
# Append some text to fastq sample name.
# Arguments:
#   Input fastq path.
#   Text to add.
# Outputs:
#   Writes modified fastq sample name to stdout.
#######################################
function append_to_sample_name() {
    sed "s/\(^.*\)_S\([[:digit:]]\)_L\([[:digit:]]\+\)_\([RI][[:digit:]]\)_001\(\.fastq\.gz\)/\1_${2}_S\2_L\3_\4_001\5/" <(basename "${1}")
}

#######################################
# Main program.
# Globals:
#   FASTQ
#   REGION
# Returns:
#   0 if completed, non-zero on error.
#######################################
function main() {
    append_to_sample_name "${FASTQ}" "${REGION}"
}

main
