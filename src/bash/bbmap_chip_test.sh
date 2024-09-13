#!/usr/bin/env bash

set -u
set -e

#######################################
# Constants
#######################################
# Array of fastq files
readonly fastqs=(reads/raw/wt_{antiflag,input}_ip1_*_001.fastq.gz)
# Path to reference genome
readonly reference_genome="references/R64-1-1/genome_I.fa"
# Path to directory where output map (BAM) files from STAR should be stored
readonly map_dir="mappings"


#######################################
# Load software
#######################################
# echo to stderr
function errcho {
    >&2 echo "${@}"
}
function bbmap {
    local image_uri="https://depot.galaxyproject.org/singularity/bbmap:39.08--h92535d8_0"
    apptainer exec "${image_uri}" bbmap.sh "${@}"
}
function samtools {
    local image_uri="https://depot.galaxyproject.org/singularity/samtools%3A1.19.2--h50ea8bc_1"
    apptainer exec "${image_uri}" samtools "${@}"
}


#######################################
# Main program.
#######################################
# create target directories
mkdir -p "${map_dir}"

# create array of fastq prefixes

## use associative array to store unique fastq path prefixes
declare -A unsorted_unique_fastq_path_prefixes
## iterate over fastqs and store path prefixes in the associative array
for fastq in "${fastqs[@]}"
do
    # extract the path prefix from the full fastq file path and add it to the associative array
    fastq_path_prefix=$(echo "${fastq}" | sed -E 's/_S1_L001_R[12]_001.fastq.gz$//g')
    unsorted_unique_fastq_path_prefixes["${fastq_path_prefix}"]=1
done

## the associative array is useful for getting unique elements
## however, it is not determinitive and is not easily indexed by number
## we need a constant ordering and easy indexing by number to run in parallel with Slurm arrays
## sort the fastq path prefixes to get a determinitive order
## store the output as a regular array which we can access by numerical index
unique_fastq_path_prefixes=($(printf '%s\n' "${!unsorted_unique_fastq_path_prefixes[@]}" | sort))

for unique_fastq_path_prefix in "${unique_fastq_path_prefixes[@]}"
do
    errcho "Fastq path prefix: ${fastq_path_prefix}"

    # get prefix for sample name
    sample_prefix=$(basename "${unique_fastq_path_prefix}")
    errcho "Sample prefix: ${sample_prefix}"

    # map reads to reference with bbmap
    bbmap \
        nodisk=true \
        ref="${reference_genome}" \
        in="${unique_fastq_path_prefix}_S1_L001_R1_001.fastq.gz" \
        in2="${unique_fastq_path_prefix}_S1_L001_R2_001.fastq.gz" \
        unpigz=true \
        fast=true \
        threads=12 \
        out=/dev/stdout \
    | samtools collate -@ 4 -O -u - \
    | samtools fixmate -@ 4 -m -u - - \
    | samtools sort -@ 4 -u - \
    | samtools markdup -@ 4 -S -d 2500 --mode s --include-fails - "${map_dir}/${sample_prefix}.bam"

    # index the BAM file
    samtools index "${map_dir}/${sample_prefix}.bam"
done

exit 0
