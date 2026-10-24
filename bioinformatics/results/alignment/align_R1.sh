#!/usr/bin/env bash
set -euo pipefail
cd /home/ravanox/codex/bioinformatics
bwa=tools/miniforge3/envs/bwa/bin/bwa
samtools=tools/samtools-env/bin/samtools
ref=data/reference/GCA_000001405.15_GRCh38_genomic.fna
out=results/alignment/571566_R1
"$bwa" mem -t 6 -K 20000000 -R '@RG\tID:571566_R1\tSM:571566' "$ref" results/alignment/571566_R1.recovered.fq 2> "${out}.bwa.log" | "$samtools" sort -@ 1 -m 256M -T "${out}.sort_tmp" -o "${out}.sorted.bam" - 2> "${out}.sort.log"
"$samtools" quickcheck -v "${out}.sorted.bam"
"$samtools" index -@ 2 "${out}.sorted.bam"
"$samtools" flagstat -@ 2 "${out}.sorted.bam" > "${out}.flagstat.txt"
"$samtools" idxstats "${out}.sorted.bam" > "${out}.idxstats.txt"
printf 'Alignment and validation completed\n'
