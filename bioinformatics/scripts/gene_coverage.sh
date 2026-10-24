#!/usr/bin/env bash
set -euo pipefail
cd /home/ravanox/codex/bioinformatics
python3 scripts/build_sfari_bed.py
for depth in 1 10 20; do
    tools/samtools-env/bin/samtools bedcov -Q 20 -G SUPPLEMENTARY -j -d "$depth" \
        data/reference/sfari.mane.coding_flanks.bed results/wgs/571566.markdup.bam \
        > "results/wgs/qc/571566.sfari_coding.min${depth}.bedcov.tsv"
done
python3 - <<'PY'
import csv
from collections import defaultdict
from pathlib import Path
qc = Path('results/wgs/qc')
genes = defaultdict(lambda: {'length':0, 'depth_sum':0, 'covered1':0, 'covered10':0, 'covered20':0})
for threshold in (1,10,20):
    for line in (qc/f'571566.sfari_coding.min{threshold}.bedcov.tsv').read_text().splitlines():
        chrom,start,end,gene,total,covered=line.split('\t')
        row=genes[gene]
        if threshold==1:
            row['length']+=int(end)-int(start)
            row['depth_sum']+=int(total)
        row[f'covered{threshold}']+=int(covered)
with (qc/'571566.sfari_gene_coverage.tsv').open('w') as handle:
    out=csv.writer(handle,delimiter='\t')
    out.writerow(['gene','coding_plus_flank_bases','mean_depth','percent_ge1','percent_ge10','percent_ge20'])
    for gene,row in sorted(genes.items()):
        n=row['length']
        out.writerow([gene,n,f"{row['depth_sum']/n:.2f}",*[f"{100*row[f'covered{x}']/n:.2f}" for x in (1,10,20)]])
PY
touch results/wgs/qc/gene_coverage.complete
