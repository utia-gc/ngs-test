#!/usr/bin/env bash

wget \
    --no-clobber \
    --directory-prefix="data/alignments" \
    --input-file="data/urls/reads_urls.txt"
