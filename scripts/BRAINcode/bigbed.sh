#!/bin/bash

conda activate /PHShome/rw552/condaenvs/ucsc

### making bigbed files
sort-bed BRAINcode.nonneuronal.TNE.hg38.bed > BRAINcode.nonneuronal.TNE.sorted.hg38.bed
bedToBigBed BRAINcode.nonneuronal.TNE.sorted.hg38.bed /PHShome/rw552/Documents/hg38.chrom.sizes BRAINcodenonneuronal.TNE.sorted.hg38.bb

