"""Union MANE coding intervals; supplement missing genes with Ensembl 115."""
import csv
from collections import defaultdict
from pathlib import Path

root = Path(__file__).resolve().parents[1]
ref = root/'data/reference'
mapping = dict(line.split() for line in (ref/'to_accessions.tsv').read_text().splitlines())
lengths = {f[0]: int(f[1]) for line in (ref/'GRCh38.fna.fai').read_text().splitlines()
           if (f := line.split('\t'))}
order = {name: i for i, name in enumerate(lengths)}
with (ref/'sfari-genes.csv').open() as handle:
    genes = {r['ensembl-id']: r['gene-symbol'] for r in csv.DictReader(handle)}
intervals = defaultdict(list)
sources = {}
def collect(path, excluded, source):
    gene = None
    with path.open() as handle:
        for line in handle:
            chrom, start, end, kind = line.rstrip('\n').split('\t')[:4]
            if kind.startswith('Chromosome_'):
                gene = None
            elif kind.startswith('Gene_'):
                gene = genes.get(kind[5:].split('.')[0])
            elif kind.startswith('Cds_') and gene and gene not in excluded and chrom in mapping:
                acc = mapping[chrom]
                intervals[(acc, gene)].append((max(0, int(start)-2), min(lengths[acc], int(end)+2)))
                sources[gene] = source
collect(ref/'mane.bed', set(), 'MANE_1.5')
if (ref/'ensembl115.bed').exists() and (ref/'ensembl115.bed').stat().st_size:
    collect(ref/'ensembl115.bed', set(sources), 'Ensembl_115_union_fallback')
rows = []
for (chrom, gene), ranges in intervals.items():
    merged = []
    for start, end in sorted(ranges):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(end, merged[-1][1])
        else:
            merged.append([start, end])
    rows.extend((chrom, start, end, gene) for start, end in merged)
rows.sort(key=lambda r: (order[r[0]], r[1], r[2], r[3]))
assert len(intervals) > 500, 'Unexpectedly few SFARI genes matched MANE coding intervals'
(ref/'sfari.mane.coding_flanks.bed').write_text(''.join('\t'.join(map(str,r))+'\n' for r in rows))
matched = {gene for chrom, gene in intervals}
(root/'results/wgs/qc/sfari_genes_without_coding_intervals.txt').write_text(
    '\n'.join(sorted(set(genes.values()) - matched))+'\n')
(root/'results/wgs/qc/sfari_gene_coverage_sources.tsv').write_text(
    'gene\ttranscript_source\n'+''.join(f'{gene}\t{source}\n' for gene,source in sorted(sources.items())))
print(f'{len(matched)} genes; {len(rows)} merged coding intervals with 2bp flanks')
