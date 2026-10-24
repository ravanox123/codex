"""Map GRCh38 resource aliases to this exact FASTA; stream a VCF to stdout."""
import argparse
import gzip
import json
import re
import sys
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('vcf')
p.add_argument('--root', default=str(Path(__file__).resolve().parents[1]))
a = p.parse_args()
root = Path(a.root)
lengths = {}
for line in (root/'data/reference/GRCh38.fna.fai').read_text().splitlines():
    fields = line.split('\t')
    lengths[fields[0]] = int(fields[1])
aliases = {name: name for name in lengths}
for line in (root/'data/reference/GRCh38_assembly_report.txt').read_text().splitlines():
    if line.startswith('#'):
        continue
    f = line.split('\t')
    if len(f) < 10 or f[4] not in lengths:
        continue
    assert int(f[8]) == lengths[f[4]], f
    for alias in (f[0], f[4], f[6], f[9]):
        if alias != 'na':
            aliases[alias] = f[4]
aliases['MT'] = aliases['M'] = aliases['chrM']
kept = skipped = 0
with gzip.open(a.vcf, 'rt') as handle:
    for line in handle:
        if line.startswith('##contig='):
            name = re.search(r'ID=([^,>]+)', line).group(1)
            length = re.search(r'(?:length|Length)=([0-9]+)', line)
            if name in aliases and length:
                assert int(length.group(1)) == lengths[aliases[name]], line
            continue
        if line.startswith('#CHROM'):
            for name, length in lengths.items():
                print(f'##contig=<ID={name},length={length}>')
        if line.startswith('#'):
            sys.stdout.write(line)
            continue
        name, remainder = line.split('\t', 1)
        if name not in aliases:
            skipped += 1
            continue
        sys.stdout.write(aliases[name] + '\t' + remainder)
        kept += 1
print(json.dumps({'resource': a.vcf, 'retained': kept,
                  'excluded_unmatched_contigs': skipped}), file=sys.stderr)
