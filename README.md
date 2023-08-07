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
