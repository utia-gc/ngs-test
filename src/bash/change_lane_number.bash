#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

#######################################
# Constants
#######################################
readonly FASTQ="${1}"
readonly LANE_NUMBER="${2}"

#######################################
# Functions
#######################################

#######################################
# Append some text to fastq sample name.
# Arguments:
#   Input fastq path.
#   Lane number, typically into a form like 'OO1'.
# Outputs:
#   Writes modified fastq sample name to stdout.
#######################################
function change_lane_number() {
    sed "s/\(^.*\)_S\([[:digit:]]\)_L\([[:digit:]]\+\)_\([RI][[:digit:]]\)_001\(\.fastq\.gz\)/\1_S\2_L${2}_\4_001\5/" <(basename "${1}")
}

#######################################
# Main program.
# Globals:
#   FASTQ
#   LANE_NUMBER
# Outputs:
#   Writes modified fastq sample name to stdout.
# Returns:
#   0 if completed, non-zero on error.
#######################################
function main() {
    change_lane_number "${FASTQ}" "${LANE_NUMBER}"
}

main
