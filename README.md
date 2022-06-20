# test-SRAlign
A central set of files for reproducible testing of SRAlign and pipelines built on SRAlign.

## Download fastq files

### Combine fastq files

```
cat \
    SRR7167175_GSM3142773_Capped_nuclear_RNA-seq_of_wild-type_embryos_replicate_1_Caenorhabditis_elegans_RNA-Seq_1.fastq.gz \
    SRR7167176_GSM3142773_Capped_nuclear_RNA-seq_of_wild-type_embryos_replicate_1_Caenorhabditis_elegans_RNA-Seq_1.fastq.gz \
    > SRR7167175-SRR7167176_GSM3142773_Capped_nuclear_RNA-seq_of_wild-type_embryos_replicate_1_Caenorhabditis_elegans_RNA-Seq_1.fastq.gz

cat \   
    SRR7167175_GSM3142773_Capped_nuclear_RNA-seq_of_wild-type_embryos_replicate_1_Caenorhabditis_elegans_RNA-Seq_2.fastq.gz \
    SRR7167176_GSM3142773_Capped_nuclear_RNA-seq_of_wild-type_embryos_replicate_1_Caenorhabditis_elegans_RNA-Seq_2.fastq.gz \
    > SRR7167175-SRR7167176_GSM3142773_Capped_nuclear_RNA-seq_of_wild-type_embryos_replicate_1_Caenorhabditis_elegans_RNA-Seq_2.fastq.gz
```
