#!/usr/bin/env bash
# Repair the failed full-transcript export before resuming sample 571566.
set -euo pipefail
cd /home/ravanox/codex/bioinformatics
out="$PWD/results/wgs"
exec 9>"$out/autism_resume.lock"
flock -n 9 || exit 1
trap 'printf "FAILED: annotation repair/resume (%s); see logs/autism_resume.log\n" "$(date -Is)" > "$out/STATUS.txt"' ERR
printf 'Rebuilding full Ensembl annotation with 5 GB heap (%s)\n' "$(date -Is)" > "$out/STATUS.txt"
tools/variants-env/bin/snpEff -Xmx5g dump -dataDir "$PWD/data/snpeff" -noLog -bed GRCh38.115 \
    > data/reference/ensembl115.bed.partial 2> "$out/logs/snpeff_full_dump_retry.log"
python3 scripts/build_wes_targets.py --annotation data/reference/ensembl115.bed.partial
mv data/reference/ensembl115.bed.partial data/reference/ensembl115.bed
printf 'Regenerating SFARI coverage with full-transcript supplement (%s)\n' "$(date -Is)" > "$out/STATUS.txt"
bash scripts/gene_coverage.sh > "$out/logs/gene_coverage_full.log" 2>&1
touch "$out/qc/gene_coverage_full.complete"
bash scripts/run_wgs.sh
