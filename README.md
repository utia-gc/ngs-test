# ngs-test diff-chip-seq

A central set of files for reproducible testing of [`utia-gc/diff-chip-seq`](https://github.com/utia-gc/diff-chip-seq).

*S. cerevisiae* is used here because of its small chromosomes.

## References

*S. cerevisiae* sample references from build R64-1-1 are stored in `data/references/R64-1-1`.

### Fetch and extract reference

Fetch references from [Illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html)

```bash
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Saccharomyces_cerevisiae_Ensembl_R64-1-1.tar.gz
tar -xzvf Saccharomyces_cerevisiae_Ensembl_R64-1-1.tar.gz
```

### Genome reference

For testing purposes, only chromosome I will be used as the reference genome:

```bash
mv Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/Chromosomes/I.fa data/references/R64-1-1/genome_I.fa
```

#### Index reference genome

```bash
singularity exec 'https://depot.galaxyproject.org/singularity/samtools%3A1.17--hd87286a_1' samtools faidx data/references/R64-1-1/genome_I.fa
```

### Transcriptome reference

#### Filter annotations GTF for features on chromosome I

```bash
awk -F '\t' '$1 == "I"' Saccharomyces_cerevisiae/Ensembl/R64-1-1/Annotation/Archives/archive-2015-07-17-14-36-40/Genes/genes.gtf > data/references/R64-1-1/annotations_I.gtf
```

#### Extract transcript sequences

Launch interactive session in Docker container with gffread

```bash
docker run -it -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/cufflinks:2.2.1--py36_2 bash
```

Inside Docker container, extract transcript sequences

```bash
gffread data/references/R64-1-1/annotations_I.gtf -g data/references/R64-1-1/genome_I.fa -w data/references/R64-1-1/transcriptome_I.fa
```

## Reads

Fetch *S. cerevisiae* RNA-seq reads and subsample them to create minimal datasets for development and testing.

### Fetch reads

```bash
nextflow run trev-f/fetch-sra-fastq -r v0.2.1 --metadata data/metadata/sra_explorer_metadata.tsv --baseDirData data/temp
```

### Sample reads

Sample reads with [`seqtk`](https://github.com/lh3/seqtk)

Create lanes for SE reads

```bash
for lane in 1 2; do
    ./src/bash/downsample_reads.bash \
        "${lane}" \
        data/temp/reads/SRR1066657_GSM1299413_WT_NR_A_Saccharomyces_cerevisiae_RNA-Seq.fastq.gz \
        50000 \
        > "data/reads/raw/SRR1066657_S3_L00${lane}_R1_001.fastq.gz"
done
```

Create lanes for PE reads

```bash
for lane in 1 2; do
    for read in 1 2; do
        ./src/bash/downsample_reads.bash \
            "${lane}" \
            "data/temp/reads/SRR6924569_GSM3073206_Saccharomyces_cerevisiae-AR_Biological_Repeat-2_Saccharomyces_cerevisiae_RNA-Seq_${read}.fastq.gz" \
            50000 \
            > "data/reads/raw/SRR6924569_S1_L00${lane}_R${read}_001.fastq.gz"
    done
done
```

### Trim reads

#### Cutadapt

```bash
./src/bash/cutadapt.sh
```

#### fastp

```bash
./src/bash/fastp.sh
```

## Align

### Align to genome with STAR

#### Build index

```bash
mkdir data/references/R64-1-1/star_index

gunzip -c data/references/R64-1-1/genome_I.fa.gz > genome.fasta
gunzip -c data/references/R64-1-1/annotations_I.gtf.gz > annotations.gtf

singularity exec docker://quay.io/biocontainers/star:2.7.10b--h9ee0642_0 \
    STAR \
    --runMode genomeGenerate \
    --genomeDir data/references/R64-1-1/star_index \
    --genomeFastaFiles genome.fasta \
    --runThreadN 3 \
    --sjdbGTFfile annotations.gtf \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbOverhang 99

rm genome.fasta annotations.gtf
```

#### Map reads

```bash
READS=( data/reads/SRR1066657_GSM1299413_WT_NR_A_Saccharomyces_cerevisiae_RNA-Seq_50000.fastq.gz data/reads/SRR1066658_GSM1299414_WT_NR_B_Saccharomyces_cerevisiae_RNA-Seq_50000.fastq.gz "data/reads/SRR6924569_GSM3073206_Saccharomyces_cerevisiae-AR_Biological_Repeat-2_Saccharomyces_cerevisiae_RNA-Seq_1_50000.fastq.gz data/reads/SRR6924569_GSM3073206_Saccharomyces_cerevisiae-AR_Biological_Repeat-2_Saccharomyces_cerevisiae_RNA-Seq_2_50000.fastq.gz" "data/reads/SRR6924589_GSM3073211_Saccharomyces_cerevisiae-AN_Biological_Repeat-1_Saccharomyces_cerevisiae_RNA-Seq_1_50000.fastq.gz data/reads/SRR6924589_GSM3073211_Saccharomyces_cerevisiae-AN_Biological_Repeat-1_Saccharomyces_cerevisiae_RNA-Seq_2_50000.fastq.gz" )

for READ in "${READS[@]}";
do
    BASE="${READ##*/}"
    BASE="${BASE%.fastq.gz}"
    echo "Align: ${BASE}"

    singularity exec docker://quay.io/biocontainers/star:2.7.10b--h9ee0642_0 \
        STAR \
            --genomeDir data/references/R64-1-1/star_index \
            --readFilesIn ${READ} \
            --runThreadN 3 \
            --outFileNamePrefix data/alignments/${BASE} \
            --readFilesCommand zcat \
            --outSAMtype BAM Unsorted
done
```

#### Sort reads by position

```bash
for BAM in data/alignments/*.bam;
do
    BASE="${BAM##*/}"
    BASE="${BAM%Aligned.out.bam}"
    echo "Sort: ${BASE}"

    singularity exec docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0 \
        samtools sort \
            -@ 3 \
            ${BAM} \
            > ${BASE}.sorted.bam
done
```

## Create large reference genome index

Total number of bases for large reference genomes may be output in exponential notation.
The pipeline cannot use this exponential notation in its sequencing depth calculations.

Produce a genome fasta index from a large reference behavior for testing purposes to fix this bug.

## Prepare ChIP-seq data

Steps for preparing ChIP-seq data can be found in the [docs](docs/chip_test_data.md)
