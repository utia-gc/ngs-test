# ngs-test

## Purpose

Provides a central set of files for reproducible testing of [`ngs`](https://github.com/utia-gc/ngs/) and pipelines built on `ngs`.

## Recommended usage

### Repository structure and name

**Create a new branch for each pipeline that is forked off of `ngs`.**
This will keep us proliferating repositories by requiring both a forked pipeline and a forked testing repo for each new type of analysis.

**Give each branch the same name as the pipeline that it is being used to test.**
For example, if a new RNA-seq pipeline called `rnaseq` is being tested, create a new branch in this repo called `rnaseq`.
