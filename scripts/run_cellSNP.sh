#!/bin/bash -e

mkdir -p cellSNP
cut -f 2 -d ',' cellranger_${sample}/outs/per_sample_outs/cellranger_${sample}/count/sample_barcodes.csv > ${sample}_barcodes.txt
cellsnp-lite -s cellranger_${sample}/outs/per_sample_outs/cellranger_${sample}/count/sample_alignments.bam -b ${sample}_barcodes.txt -O cellSNP/${sample} -R genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf -p 8 --minMAF 0.1 --minCOUNT 20
