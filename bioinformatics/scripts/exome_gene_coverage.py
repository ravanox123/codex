"""Measure per-gene primary MANE CDS+2bp coverage with one sequential bedcov pass."""
import csv,collections,subprocess,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/wgs/exome_review'
aliases=dict(l.split() for l in (ROOT/'data/reference/to_accessions.tsv').read_text().splitlines())
lengths={f[0]:int(f[1]) for l in (ROOT/'data/reference/GRCh38.fna.fai').read_text().splitlines() if (f:=l.split())};order={c:i for i,c in enumerate(lengths)}
intervals=collections.defaultdict(list)
with (ROOT/'data/reference/mane.features.tsv').open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['type']=='Cds' and r['chr'] in aliases and r['geneName']:
   c=aliases[r['chr']];intervals[(c,r['geneName'])].append((max(0,int(r['start'])-2),min(lengths[c],int(r['end'])+2)))
rows=[]
for (chrom,gene),iv in intervals.items():
 merged=[]
 for a,b in sorted(iv):
  if merged and a<=merged[-1][1]:merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
  else:merged.append((a,b))
 rows.extend((chrom,a,b,gene) for a,b in merged)
rows.sort(key=lambda r:(order[r[0]],r[1],r[2],r[3]));bed=OUT/'mane_all_coding_flanks.bed';bed.write_text(''.join('\t'.join(map(str,r))+'\n' for r in rows))
print('Coverage intervals:',len(rows),'gene/chromosome groups:',len(intervals),flush=True)
raw=OUT/'mane_all_coding.min20.bedcov.tsv'
with raw.with_suffix('.partial').open('w') as f:
 subprocess.run([str(ROOT/'tools/samtools-env/bin/samtools'),'bedcov','-Q','20','-G','SUPPLEMENTARY','-j','-d','20',str(bed),str(ROOT/'results/wgs/571566.markdup.bam')],stdout=f,check=True)
raw.with_suffix('.partial').replace(raw)
genes=collections.defaultdict(lambda:[0,0,0])
for l in raw.read_text().splitlines():
 c,a,b,g,total,covered=l.split('\t');v=genes[g];v[0]+=int(b)-int(a);v[1]+=int(total);v[2]+=int(covered)
output=[]
for g,(n,total,covered) in sorted(genes.items()):output.append(dict(gene=g,coding_plus_flank_bases=n,mean_depth=f'{total/n:.2f}',percent_ge20=f'{100*covered/n:.2f}'))
for name,rs in [('all_gene_coverage.tsv',output),('all_gene_coverage_gaps.tsv',[r for r in output if float(r['percent_ge20'])<90])]:
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,output[0].keys(),delimiter='\t');w.writeheader();w.writerows(rs)
(OUT/'coverage_summary.json').write_text(json.dumps({'measured_genes':len(output),'genes_ge90percent_at20x':sum(float(r['percent_ge20'])>=90 for r in output),'genes_below90percent_at20x':sum(float(r['percent_ge20'])<90 for r in output),'genes_with_zero_at20x':sum(float(r['percent_ge20'])==0 for r in output)},indent=2)+'\n')
print((OUT/'coverage_summary.json').read_text(),flush=True)
