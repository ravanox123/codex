#!/usr/bin/env bash
set -euo pipefail
cd /home/ravanox/codex/bioinformatics
bcf="$PWD/tools/variants-env/bin/bcftools"
ref="$PWD/data/reference/GRCh38.fna"
for resource in 1000G_phase1.snps.high_confidence.hg38 Mills_and_1000G_gold_standard.indels.hg38 clinvar.GRCh38; do
    output="data/reference/${resource}.matched.vcf.gz"
    if [[ ! -f "${output}.tbi" ]]; then
        python3 scripts/remap_resource_vcf.py "data/reference/${resource}.vcf.gz" \
            | "$bcf" norm -f "$ref" -c x -Ou \
            | "$bcf" sort -m 384M -T results/wgs/tmp/resource-sort -Oz -o "${output}.partial.gz"
        "$bcf" index -t "${output}.partial.gz"
        mv "${output}.partial.gz" "$output"
        mv "${output}.partial.gz.tbi" "${output}.tbi"
    fi
done
date -Is > results/wgs/resources.complete
