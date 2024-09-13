# ChIP-seq Test Data

This document details how ChIP-seq test files were produced.

## Fetch raw reads from SRA

Raw reads for 2 replicates each of WT ChIP (anti-FLAG) and WT Input ChIP-seq libraries from BioProject [PRJNA1074743](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1074743) were fetched using my custom fastq fetching pipeline.

```bash
nextflow run trev-f/fetch-sra-fastq \
    -r v0.2.1 \
    --metadata metadata/sra_explorer_metadata_chip.txt
```

## Map reads to reference

Map raw reads to chromosome I of the *S. cerevisiae* reference R64-1-1.
Use BBMap for its speed and ease of use.

```bash
bash src/bash/bbmap_chip.sh
```

## Sample reads and convert to fastq format

```bash
bash src/bash/sample_chip_reads.sh
```

## Produce test BAM files

```bash
bash src/bash/bbmap_chip_test.sh
```
