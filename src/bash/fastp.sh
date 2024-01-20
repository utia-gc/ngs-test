#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

#######################################
# Aliases
#######################################
function fastp() {
    singularity exec https://depot.galaxyproject.org/singularity/fastp%3A0.23.4--hadf994f_2 fastp "${@}"
}
export -f fastp

#######################################
# Constants
#######################################
# directory for output reads and reports
readonly OUT_DIR="data/reads/trimmed/fastp"

#######################################
# Functions
#######################################

#######################################
# Run fastp for single end reads.
# Globals:
#   OUT_DIR
# Arguments:
#   Read 1 fastq path.
#   Stem name.
#######################################
function fastp_se() {
    local r1="${1}"
    local stem_name="${2}"

    fastp \
        --in1 "${r1}" \
        --out1 "${OUT_DIR}/${stem_name}"_trimmed_R1.fastq.gz \
        --json "${OUT_DIR}/${stem_name}"_fastp.json \
        --html ""
}

#######################################
# Run fastp for paired end reads.
# Globals:
#   OUT_DIR
# Arguments:
#   Read 1 fastq path.
#   Read 2 fastq path.
#   Stem name.
#######################################
function fastp_pe() {
    local r1="${1}"
    local r2="${2}"
    local stem_name="${3}"

    fastp \
        --in1 "${r1}" \
        --in2 "${r2}" \
        --out1 "${OUT_DIR}/${stem_name}"_trimmed_R1.fastq.gz \
        --out2 "${OUT_DIR}/${stem_name}"_trimmed_R2.fastq.gz \
        --json "${OUT_DIR}/${stem_name}"_fastp.json \
        --html ""
}

#######################################
# Generate stem name from an R1 file
# Arguments:
#   Read 1 fastq path.
#######################################
function generate_stem_name() {
    local r1="${1}"

    basename "${r1}" _R1_001.fastq.gz \
        | sed 's/_S[0-9]\+_/_/'
}

#######################################
# Generate the path to an R2 file from an R1 file
# Arguments:
#   Read 1 fastq path.
#######################################
function generate_r2_path() {
    local r1="${1}"

    echo "${r1}" \
        | sed 's/_R1_/_R2_/'
}

#######################################
# Main program.
# Globals:
#   OUT_DIR
#######################################
function main() {
    mkdir -p "${OUT_DIR}"

    for r1 in data/reads/raw/SRR*_S*_L00?_R1_001.fastq.gz; do
        echo "Parsing reads associated with: ${r1}" >&2
        
        r2=$(generate_r2_path "${r1}")
        stem_name=$(generate_stem_name "${r1}")
        
        if [ -e "${r2}" ]; then
            echo "${stem_name} has paired-end reads." >&2
            echo "Run fastp in PE mode. . ." >&2
            fastp_pe "${r1}" "${r2}" "${stem_name}"
        else
            echo "${stem_name} has single-end reads." >&2
            echo "Run fastp in SE mode. . ." >&2
            fastp_se "${r1}" "${stem_name}"
        fi
        
        echo >&2
    done
}

main
