#!/usr/bin/env bash

wget \
    --no-clobber \
    --directory-prefix="data/temp" \
    --input-file="data/urls/reads_urls.txt"
