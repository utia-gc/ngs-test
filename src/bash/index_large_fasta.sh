#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

#######################################
# Aliases
#######################################
function samtools {
    singularity exec https://depot.galaxyproject.org/singularity/samtools%3A1.17--hd87286a_1 samtools "${@}"
}
export -f samtools

#######################################
# Constants
#######################################
# directory for output refs
readonly OUT_DIR="data/references"
# url for genome fasta
readonly GENOME_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz"

#######################################
# Functions
#######################################

#######################################
# Index fasta
# Arguments:
#   Path to fasta file
#######################################
function index_fasta() {
    local fasta="${1}"

    echo "Index fasta: ${fasta}" >&2
    samtools faidx "${fasta}"
}

#######################################
# Fetch the genome file
# Globals:
#   GENOME_URL
#   OUT_DIR
#######################################
function fetch_genome() {
    echo "Fetching genome: ${GENOME_URL}" >&2

    wget \
        --no-clobber \
        --directory-prefix=${OUT_DIR} \
        "${GENOME_URL}"
}

#######################################
# Generate the path to an R2 file from an R1 file
# Globals:
#   GENOME_URL
#   OUT_DIR
#######################################
function generate_genome_path() {
    local base=$(basename "${GENOME_URL}")

    echo "${OUT_DIR}/${base}"
}

#######################################
# Main program.
# Globals:
#   OUT_DIR
#######################################
function main() {
    mkdir -p "${OUT_DIR}"

    fetch_genome

    genome_path=$(generate_genome_path)
    echo "Path to genome: ${genome_path}" >&2

    echo "Unzip fasta" >&2
    unpigz "${genome_path}"

    # use parameter expansion to remove the .gz extension since the file has been decompressed
    unzipped_genome_path="${genome_path%.gz}"
    index_fasta "${unzipped_genome_path}"

    # cleanup the fasta file if the index is made
    if [ -e "${unzipped_genome_path}.fai" ]; then
        echo "Remove unzipped fasta: ${unzipped_genome_path}" >&2
        rm "${unzipped_genome_path}"
    fi
}

main
