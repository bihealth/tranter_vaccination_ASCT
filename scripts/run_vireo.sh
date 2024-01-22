#!/bin/bash

unset DISPLAY
[ ! -f cellSNP/${sample}/cellSNP.base.vcf.gz ] && gzip cellSNP/${sample}/cellSNP.base.vcf
[ ! -f cellSNP/${sample}/cellSNP.samples.tsv.gz ] && gzip cellSNP/${sample}/cellSNP.samples.tsv
vireo -c cellSNP/${sample} -N ${N} -o vireo/${sample}
GTbarcode -i vireo/${sample}/GT_donors.vireo.vcf.gz -o vireo/${sample}/GT_barcodes.tsv --randSeed 1
