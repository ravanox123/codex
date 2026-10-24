"""Build approximate WES calling targets from primary-chromosome exons."""
import argparse
from collections import defaultdict
from pathlib import Path

root = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument('--annotation', type=Path, default=root/'data/reference/ensembl115.bed')
p.add_argument('--output', type=Path, default=root/'data/reference/wes.exons_padded100.bed')
a = p.parse_args()
ref = root/'data/reference'
aliases = dict(line.split() for line in (ref/'to_accessions.tsv').read_text().splitlines())
lengths = {f[0]: int(f[1]) for line in (ref/'GRCh38.fna.fai').read_text().splitlines()
           if (f := line.split('\t'))}
ranges = defaultdict(list)
with a.annotation.open() as handle:
    for line in handle:
        chrom, start, end, feature = line.rstrip('\n').split('\t')[:4]
        if chrom in aliases and feature.startswith('Exon_'):
            acc = aliases[chrom]
            ranges[acc].append((max(0, int(start)-100), min(lengths[acc], int(end)+100)))
assert sum(map(len, ranges.values())) > 10000, 'Incomplete exon annotation'
rows = []
for chrom in lengths:
    merged = []
    for start, end in sorted(ranges.get(chrom, [])):
        assert 0 <= start < end <= lengths[chrom]
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    rows.extend((chrom, start, end) for start, end in merged)
assert len({r[0] for r in rows}) == 24, 'Expected all nuclear primary chromosomes'
temporary = a.output.with_name(a.output.name+'.partial')
temporary.write_text(''.join(f'{chrom}\t{start}\t{end}\n' for chrom,start,end in rows))
temporary.replace(a.output)
print(f'Approximate WES targets: {len(rows)} merged exon intervals; '
      f'{sum(end-start for _,start,end in rows):,} bp including 100bp flanks. '
      'These are annotation-derived targets, not the capture-kit BED.')
