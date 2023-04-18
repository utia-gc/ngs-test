# test-SRAlign - RNA-seq

A central set of files for reproducible testing of RNA-seq pipelines built on SRAlign.

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
nextflow run trev-f/fetch-sra-fastq -r v0.2.1 --metadata data/metadata/sra_explorer_metadata.tsv
```

### Sample reads

Launch interactive session in Docker container with seqtk.

```bash
docker run -it -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/seqtk:1.3--h7132678_4 bash
```

Inside Docker container, sample reads.

```bash
for FQ in data/reads/raw/*.fastq.gz;
do
    echo "Sampling: ${FQ}"
    BASE="${FQ##*/}"
    BASE="${BASE%.fastq.gz}"
    seqtk sample -s100 ${FQ} 50000 | gzip -c > data/reads/${BASE}_50000.fastq.gz
done
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
