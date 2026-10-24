"""Independent CIGAR/base support summaries at the seven ClinVar assertion SNVs."""
from pathlib import Path
import json,subprocess,re,collections,statistics
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/wgs/exome_review';aliases=dict(l.split() for l in (ROOT/'data/reference/to_accessions.tsv').read_text().splitlines())
rows=[d['candidate'] for d in json.loads((OUT/'details.json').read_text()) if d['candidate']['review_tier']=='1_explicit_ClinVar_P_LP_review'];output=[]
for r in rows:
 pos=int(r['position']);chrom=aliases[r['chromosome']];sam=subprocess.check_output([str(ROOT/'tools/samtools-env/bin/samtools'),'view','-F','3844',str(ROOT/'results/wgs/571566.recal.bam'),f'{chrom}:{pos}-{pos}'],text=True)
 counts=collections.Counter();altmqs=[];altnms=[];refnms=[];starts=set()
 for line in sam.splitlines():
  c=line.split('\t');rp=int(c[3]);qp=0;base=None;quality=0
  for n,op in re.findall(r'(\d+)([MIDNSHP=X])',c[5]):
   n=int(n)
   if op in 'M=X':
    if rp<=pos<rp+n:
     off=qp+pos-rp;base=c[9][off];quality=ord(c[10][off])-33
    rp+=n;qp+=n
   elif op in 'DN':rp+=n
   elif op in 'IS':qp+=n
  if base is None or quality<20:continue
  kind='alt' if base==r['alt'] else 'ref' if base==r['ref'] else 'other';mq=int(c[4]);nm=next((int(t[5:]) for t in c[11:] if t.startswith('NM:i:')),None)
  counts[kind+'_BQ20_anyMQ']+=1
  if mq>=20:
   counts[kind+'_BQ20_MQ20']+=1;counts[kind+('_reverse' if int(c[1])&16 else '_forward')]+=1
   if kind=='alt':
    altmqs.append(mq);starts.add((c[3],c[5],bool(int(c[1])&16)))
    if nm is not None:altnms.append(nm)
   if kind=='ref' and nm is not None:refnms.append(nm)
  if mq>=40:counts[kind+'_BQ20_MQ40']+=1
  if mq>=60:counts[kind+'_BQ20_MQ60']+=1
 x={'gene':r['gene'],'chromosome':r['chromosome'],'position':pos,'ref':r['ref'],'alt':r['alt'],'caller_GT':r['genotype'],'caller_AD':r['allele_depths'],**dict(counts),'alt_mean_mapq':statistics.mean(altmqs) if altmqs else None,'alt_median_NM':statistics.median(altnms) if altnms else None,'ref_median_NM':statistics.median(refnms) if refnms else None,'alt_distinct_start_cigar_strand':len(starts)}
 output.append(x);print(x,flush=True)
(OUT/'assertion_SNV_read_support.json').write_text(json.dumps(output,indent=2)+'\n')
