# ngs-test

## Purpose

Provides a central set of files for reproducible testing of [`ngs`](https://github.com/utia-gc/ngs/) and pipelines built on `ngs`.

## Recommended usage

### Repository structure and name

**Create a new branch for each pipeline that is forked off of `ngs`.**
This will keep us proliferating repositories by requiring both a forked pipeline and a forked testing repo for each new type of analysis.

**Give each branch the same name as the pipeline that it is being used to test.**
For example, if a new RNA-seq pipeline called `rnaseq` is being tested, create a new branch in this repo called `rnaseq`.

## Test files

### Good testing practices

**Use test files that are similar to the real data as possible.**
This seems obvious, but it's worth emphasizing.
This means that if the data you commonly work with typically comes straight from the sequencing core, then use names similar to what they use.
If your data typically comes from a specific kit version, use sample data that comes from that kit.
The closer the tests data resembles the real data your pipeline will work with, the fewer bugs you're likely to encounter.

### Getting and processing data

**Use downsized data wherever possible.**
Ideally, testing should be quick and straightforward to allow for fast iteration during development.
Using full datasets which can take several minutes if not hours to finish processing is not amenable to this goal of quick iteration.
Additionally, standard GitHub repos do not allow for the storage of files in excess of 100 MB.
Unless for some rare case it is genuinely impossible to use downsized data, we cannot recommend downsizing enough.

**Downsize data by sampling from a single, small chromosome.**
Our preferred way to downsize data is to sample data only from a single chromosome.
Below are some strategies for how we accomplish this.

* Reference genomes: Download reference sequence for single chromosome. Most reference genome repositories will host download options for single chromosome fasta files.
* Reference annotations: Both GFF and GTF annotation formats are tab-delimited file formats with the chromosome name in the first field of each line. It is trivially easy to filter these files for features on a specific chromosome using an `awk` one liner:

```bash
awk -F '\t' '$1 == "<chromosome name>"' annotations.gtf > annotations_<chromsome name>.gtf
```
