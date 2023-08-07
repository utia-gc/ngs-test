# ngs-test - scrnaseq

## Purpose

Provides a central set of files for reproducible testing of `scrnaseq`.

## Data source

### Fetch reads data

Read data was downloaded from 10x Genomics datasets page for dataset [5k Human PBMCs, 3' v3.1, Chromium Controller](https://www.10xgenomics.com/resources/datasets/5k-human-pbmcs-3-v3-1-chromium-controller-3-1-standard).
Reads were downloaded in genome-aligned BAM format (along with BAM index) to allow for easy extraction of reads that only map to a single chromosome during downsizing and downsampling.

```bash
bash ./src/bash/fetch_reads.bash
```
