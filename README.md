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

### Sample reads

Sampled to 100,000 reads each end

```
seqtk sample -s100 SRR7167175-SRR7167176_GSM3142773_Capped_nuclear_RNA-seq_of_wild-type_embryos_replicate_1_Caenorhabditis_elegans_RNA-Seq_1.fastq.gz 100000 | gzip -c > SRR7167175-SRR7167176_GSM3142773_Capped_nuclear_RNA-seq_of_wild-type_embryos_replicate_1_Caenorhabditis_elegans_RNA-Seq__skS-100000_1.fastq.gz

seqtk sample -s100 SRR7167175-SRR7167176_GSM3142773_Capped_nuclear_RNA-seq_of_wild-type_embryos_replicate_1_Caenorhabditis_elegans_RNA-Seq_2.fastq.gz 100000 | gzip -c > SRR7167175-SRR7167176_GSM3142773_Capped_nuclear_RNA-seq_of_wild-type_embryos_replicate_1_Caenorhabditis_elegans_RNA-Seq__skS-100000_2.fastq.gz

seqtk sample -s100 SRR7167177_GSM3142774_Capped_nuclear_RNA-seq_of_wild-type_embryos_replicate_2_Caenorhabditis_elegans_RNA-Seq_1.fastq.gz 100000 | gzip -c > SRR7167177_GSM3142774_Capped_nuclear_RNA-seq_of_wild-type_embryos_replicate_2_Caenorhabditis_elegans_RNA-Seq__skS-100000_1.fastq.gz

seqtk sample -s100 SRR7167177_GSM3142774_Capped_nuclear_RNA-seq_of_wild-type_embryos_replicate_2_Caenorhabditis_elegans_RNA-Seq_2.fastq.gz 100000 | gzip -c > SRR7167177_GSM3142774_Capped_nuclear_RNA-seq_of_wild-type_embryos_replicate_2_Caenorhabditis_elegans_RNA-Seq__skS-100000_2.fastq.gz
```
