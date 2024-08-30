#!/usr/bin/env bash

set -u
set -e

#######################################
# Constants
#######################################
# Array of SAM files of mapped reads
readonly sams=(data/mappings/*_mapped.sam)
# Target number of records (read pairs)
readonly target_n_records=100000
# Subsample reads seed. should be between 0-9
readonly seed="1"
# Directory to write output reads into
readonly reads_dir="reads/raw"


#######################################
# Load software
#######################################
# echo to stderr
function errcho {
    >&2 echo "${@}"
}
function samtools {
    local image_uri="https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1"
    apptainer exec "${image_uri}" samtools "${@}"
}


#######################################
# Main program.
#######################################
# create target directories
for sam in "${sams[@]}"
do
    errcho "SAM file path: ${sam}"

    # construct a simpler sample name
    sample_prefix=$(basename "${sam}" .sam | sed -E 's/^SRR[0-9]+_GSM[0-9]+_//g;s/_Saccharomyces_cerevisiae_Kluyveromyces_lactis_ChIP-Seq_mapped$//g' | tr '[:upper:]' '[:lower:]')
    errcho "Simple sample prefix: ${sample_prefix}"

    # compute the subsampling proportion
    existing_n_records=$(samtools view -c "${sam}")
    errcho "Number records in SAM: ${existing_n_records}"
    subsample_proportion=$(python3 -c "print(${target_n_records} / ${existing_n_records})")
    errcho "Subsampling proportion: ${subsample_proportion}"

    # sample the bams and write to fastqs
    samtools view \
        -uh \
        --subsample "${subsample_proportion}" \
        --subsample-seed "${seed}" \
        "${sam}" | \
    samtools fastq \
        -N \
        -1 "${reads_dir}/${sample_prefix}_S1_L00${seed}_R1_001.fastq.gz" \
        -2 "${reads_dir}/${sample_prefix}_S1_L00${seed}_R2_001.fastq.gz" \
        -
done
