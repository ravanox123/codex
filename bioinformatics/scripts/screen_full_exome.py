"""Phenotype-free audit of all called variants; preserve unreviewed calls explicitly."""
import csv,gzip,json,re,collections,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/wgs/exome_review';OUT.mkdir(exist_ok=True)
names=dict(l.split() for l in (ROOT/'data/reference/to_standard.tsv').read_text().splitlines());mane=set()
with (ROOT/'data/reference/mane.features.tsv').open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['type']=='Transcript':mane.add(r['id'].split('.')[0])
counts=collections.Counter();rows=[];details=[]
with gzip.open(ROOT/'results/wgs/571566.annotated.vcf.gz','rt') as f:
 for l in f:
  if l.startswith('#'):continue
  v=l.rstrip().split('\t');counts['all_variants']+=1;i=dict(x.split('=',1) if '=' in x else (x,True) for x in v[7].split(';'));fmt=dict(zip(v[8].split(':'),v[9].split(':')))
  if not any(x not in ('0','.') for x in fmt.get('GT','.').replace('|','/').split('/')):continue
  anns=[a.split('|') for a in i.get('ANN','').split(',') if len(a.split('|'))>10 and a.split('|')[0]==v[4]]
  hi=[a for a in anns if a[2]=='HIGH'];mh=[a for a in hi if a[6].split('.')[0] in mane];sig=i.get('CLNSIG','');tok=set(re.split(r'[/|,]',sig));p=bool(tok&{'Pathogenic','Likely_pathogenic'});conflict='Conflicting' in sig
  passed=v[6]=='PASS';counts['passing']+=passed;counts['passing_HIGH_any']+=passed and bool(hi);counts['passing_HIGH_MANE']+=passed and bool(mh);counts['passing_explicit_P_LP']+=passed and p
  if not(p or hi):continue
  selected=mh or [a for a in anns if a[6].split('.')[0] in mane] or hi or anns
  if not selected:continue
  a=selected[0];chrom=names.get(v[0],v[0]);sex=chrom in ('X','Y')
  try:af=max(float(x) for x in i.get('POP_AF','.').split(','))
  except ValueError:af=None
  ad=fmt.get('AD','.');fraction=None
  try:d=list(map(int,ad.split(',')));fraction=d[1]/sum(d) if sum(d) else None
  except (ValueError,IndexError):pass
  het='0' in fmt.get('GT','').replace('|','/').split('/') and '1' in fmt.get('GT','').replace('|','/').split('/')
  notes=[]
  if not passed:notes.append('fails original caller/genotype filters')
  if sex:notes.append('sex-chromosome genotype/ploidy unresolved')
  if af is None:notes.append('population frequency missing; rarity unproven')
  elif af>=.01:notes.append('common in limited frequency resource')
  if het and (fraction is None or not .2<=fraction<=.8):notes.append('heterozygous allele imbalance')
  if hi and not mh:notes.append('HIGH effect absent from matched MANE annotations')
  if tok&{'Benign','Likely_benign'}:notes.append('benign/likely-benign ClinVar assertion')
  if 'WARNING' in i.get('ANN','') or 'ERROR' in i.get('ANN',''):notes.append('annotation warning/error')
  if not passed:tier='5_filtered_call'
  elif sex:tier='4_sex_ploidy_unresolved'
  elif p:tier='1_explicit_ClinVar_P_LP_review'
  elif af is not None and af>=.01 or tok&{'Benign','Likely_benign'}:tier='4_known_common_or_benign'
  elif mh:tier='2_MANE_HIGH_unreviewed'
  else:tier='3_alternative_HIGH_unreviewed'
  r=dict(chromosome=chrom,position=v[1],ref=v[3],alt=v[4],gene=a[3],effect=a[1],predicted_impact=a[2],transcript=a[6],hgvs_c=a[9],hgvs_p=a[10],genotype=fmt.get('GT','.'),depth=fmt.get('DP','.'),genotype_quality=fmt.get('GQ','.'),allele_depths=ad,alt_fraction='.' if fraction is None else f'{fraction:.3f}',filter=v[6],clinvar_significance=sig or '.',clinvar_condition=i.get('CLNDN','.'),clinvar_review=i.get('CLNREVSTAT','.'),clinvar_allele_id=i.get('ALLELEID','.'),population_af=i.get('POP_AF','.'),review_tier=tier,review_notes='; '.join(notes),all_HIGH_genes=';'.join(sorted({a[3] for a in hi})),mane_HIGH_genes=';'.join(sorted({a[3] for a in mh})))
  rows.append(r);details.append({'candidate':r,'info':{k:x for k,x in i.items() if k!='ANN'},'annotations':anns})
rows.sort(key=lambda r:(r['review_tier'],r['gene'],int(r['position'])))
# Nearby indels can jointly restore reading frame; grouping alone does not establish phase.
bygene=collections.defaultdict(list)
for r in rows:
 if len(r['ref'])!=len(r['alt']):bygene[(r['chromosome'],r['gene'])].append(r)
for group in bygene.values():
 group.sort(key=lambda r:int(r['position']))
 for r in group:
  neighbors=[n for n in group if n is not r and abs(int(n['position'])-int(r['position']))<=50]
  if neighbors:r['review_notes']+='; nearby same-gene indel(s); joint haplotype/frame requires review'
fields=list(rows[0])
def write(name,items):
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(items)
write('all_HIGH_or_pathogenic_assertion.tsv',rows)
write('explicit_pathogenic_assertions.tsv',[r for r in rows if r['review_tier']=='1_explicit_ClinVar_P_LP_review'])
write('mane_HIGH_candidates.tsv',[r for r in rows if r['review_tier']=='2_MANE_HIGH_unreviewed'])
# Retain all variants in genes of key clinical-assertion leads for second-allele/nearby-call review.
keygenes={'FDXR','WNT10A','VWF','DNAH14','CNN2','SMARCA5','CNTN5','PRKAA2','DOLK','TRMT5','PSPH','SUZ12','RBMX'}
# Details for all high/assertion calls support additional targeted read checks.
(OUT/'details.json').write_text(json.dumps(details,indent=2)+'\n')
keydetails=[d for d in details if d['candidate']['gene'] in keygenes and d['candidate']['filter']=='PASS' and d['candidate']['review_tier']!='4_known_common_or_benign']
(OUT/'key_details.json').write_text(json.dumps(keydetails,indent=2)+'\n')
counts['all_HIGH_or_P_LP_records']=len(rows);counts['unreviewed_MANE_HIGH_candidates']=sum(r['review_tier']=='2_MANE_HIGH_unreviewed' for r in rows)
(OUT/'summary.json').write_text(json.dumps(counts,indent=2)+'\n')
print(json.dumps(counts,indent=2));print('Key read checks:',len(keydetails))
