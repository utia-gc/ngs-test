#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

#######################################
# Aliases
#######################################
function cutadapt() {
    singularity exec https://depot.galaxyproject.org/singularity/cutadapt%3A4.4--py39hf95cd2a_1 cutadapt "${@}"
}
export -f cutadapt

#######################################
# Constants
#######################################
# directory for output reads and reports
readonly OUT_DIR="data/reads/trimmed/cutadapt"
# adapter for read 1
readonly R1_ADAPTER="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
# adapter for read 2
readonly R2_ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
# minimum length for passing reads
readonly MIN_LEN="30"

#######################################
# Functions
#######################################

#######################################
# Run cutadapt for single end reads.
# Globals:
#   OUT_DIR
#   R1_ADAPTER
#   MIN_LEN
# Arguments:
#   Read 1 fastq path.
#   Stem name.
#######################################
function cutadapt_se() {
    local r1="${1}"
    local stem_name="${2}"

    cutadapt \
        -a "${R1_ADAPTER}" \
        -m "${MIN_LEN}" \
        -o "${OUT_DIR}/${stem_name}"_trimmed_R1.fastq.gz \
        "${r1}" \
        > "${OUT_DIR}/${stem_name}"_cutadapt-log.txt
}

#######################################
# Run cutadapt for paired end reads.
# Globals:
#   OUT_DIR
#   R1_ADAPTER
#   R2_ADAPTER
#   MIN_LEN
# Arguments:
#   Read 1 fastq path.
#   Read 2 fastq path.
#   Stem name.
#######################################
function cutadapt_pe() {
    local r1="${1}"
    local r2="${2}"
    local stem_name="${3}"

    cutadapt \
        -a "${R1_ADAPTER}" \
        -A "${R2_ADAPTER}" \
        -m "${MIN_LEN}" \
        -o "${OUT_DIR}/${stem_name}"_trimmed_R1.fastq.gz \
        -p "${OUT_DIR}/${stem_name}"_trimmed_R2.fastq.gz \
        "${r1}" "${r2}" \
        > "${OUT_DIR}/${stem_name}"_cutadapt-log.txt
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
            echo "Run cutadapt in PE mode. . ." >&2
            cutadapt_pe "${r1}" "${r2}" "${stem_name}"
        else
            echo "${stem_name} has single-end reads." >&2
            echo "Run cutadapt in SE mode. . ." >&2
            cutadapt_se "${r1}" "${stem_name}"
        fi
        
        echo >&2
    done
}

main
