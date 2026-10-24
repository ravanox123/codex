"""Conservative local haplotype support check; not a clinical validation or IGV review."""
import collections,csv,json,re,subprocess,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--details',type=Path,default=ROOT/'results/wgs/candidate_review/details.json');parser.add_argument('--output',type=Path,default=ROOT/'results/wgs/candidate_review/read_support.json');args=parser.parse_args()
OUT=args.output.parent
SAM=str(ROOT/'tools/samtools-env/bin/samtools')
BAM=str(ROOT/'results/wgs/571566.recal.bam')
REF=str(ROOT/'data/reference/GRCh38.fna')
aliases=dict(l.split() for l in (ROOT/'data/reference/to_accessions.tsv').read_text().splitlines())
results=[]
for d in json.loads(args.details.read_text()):
 r=d['candidate'];pos=int(r['position']);chrom=aliases[r['chromosome']];start=pos-25;end=pos+len(r['ref'])-1+25
 fasta=subprocess.check_output([SAM,'faidx',REF,f'{chrom}:{start}-{end}'],text=True)
 ref=''.join(fasta.splitlines()[1:]).upper();assert ref[25:25+len(r['ref'])]==r['ref']
 alt=ref[:25]+r['alt']+ref[25+len(r['ref']):]
 sam=subprocess.check_output([SAM,'view','-q','20','-F','3844',BAM,f'{chrom}:{start}-{end}'],text=True)
 counts=collections.Counter();strands=collections.Counter();starts=collections.defaultdict(set);mqs=collections.defaultdict(list);other=collections.Counter()
 for line in sam.splitlines():
  c=line.split('\t');flag=int(c[1]);left=int(c[3]);rp=left;qp=0;parts=[];quals=[];covered=set();bad=False
  for length,op in re.findall(r'(\d+)([MIDNSHP=X])',c[5]):
   n=int(length)
   if op in 'M=X':
    for k in range(n):
     if start<=rp+k<=end:
      parts.append(c[9][qp+k]);quals.append(ord(c[10][qp+k])-33 if c[10]!='*' else 0);covered.add(rp+k)
    rp+=n;qp+=n
   elif op=='I':
    if start<=rp-1<end:
     parts.append(c[9][qp:qp+n]);quals.extend(ord(x)-33 for x in c[10][qp:qp+n])
    qp+=n
   elif op in 'DN':
    if op=='N' and rp<=end and rp+n>start:bad=True
    covered.update(range(max(start,rp),min(end+1,rp+n)));rp+=n
   elif op=='S':qp+=n
  if bad or len(covered)!=(end-start+1) or not quals or min(quals)<20:
   counts['excluded_partial_or_low_quality']+=1;continue
  seq=''.join(parts).upper();kind='ref' if seq==ref else 'alt' if seq==alt else 'other'
  if kind=='other':
   # Allow up to two flank substitutions, preserving the exact focal allele.
   # This accommodates nearby SNPs without treating other repeat lengths as support.
   for label, expected in (('ref',ref),('alt',alt)):
    if len(seq)==len(expected) and seq[25:-25]==expected[25:-25]:
     flank_diffs=sum(a!=b for a,b in zip(seq[:25]+seq[-25:],expected[:25]+expected[-25:]))
     if 0<flank_diffs<=2:
      counts[label+'_with_flank_substitutions']+=1;kind=label;break
  counts[kind]+=1;strands[kind+('_reverse' if flag&16 else '_forward')]+=1;starts[kind].add((left,c[5],bool(flag&16)));mqs[kind].append(int(c[4]))
  if kind=='other':other[seq]+=1
 result={'gene':r['gene'],'chromosome':r['chromosome'],'position':pos,'reference_allele':r['ref'],'alternate_allele':r['alt'],**dict(counts),**dict(strands),'alt_distinct_start_cigar_strand':len(starts['alt']),'alt_mean_mapq':round(sum(mqs['alt'])/len(mqs['alt']),1) if mqs['alt'] else None,'other_haplotypes':dict(other)}
 results.append(result);print(r['gene'],dict(counts),dict(strands),flush=True)
args.output.write_text(json.dumps(results,indent=2)+'\n')
