"""Assemble the documented manual triage of the original 18 research candidates."""
import csv,json,re,hashlib
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/wgs/candidate_review'
IDS=dict(USH2A=48427,EXOC6B=1601386,DST=809956,ENPP1=13589,DLX6=1385415,STXBP1=365052,NOTCH1=477890,POLR3A=368925,PTEN=440216,SYCE1=2672423,CTR9=1091039,VDR=1252227,BRCA2=51231,MYO1E=983002,ZFHX3=808107,POLR2A=1686168,SETBP1=326725,KATNAL2=3060277)
# These are review dispositions, not independent ACMG classifications.
NOTES={
'USH2A':('Lower priority; benign-leaning conflict','Most contributing submissions are benign/likely benign, with one uncertain. ClinVar lists gnomAD frequencies around 0.42–0.68%, varying by dataset. This heterozygous change does not establish a recessive USH2A diagnosis.','Reconsider only with relevant hearing/retinal findings and a full inheritance-aware analysis.'),
'EXOC6B':('Deprioritize; benign transcript-dependent flag','The MANE transcript places the duplication in an intron (c.1980+15_1980+16dup), rather than the alternative-transcript frameshift selected by the screen. ClinVar is benign; population datasets listed there place the allele near 99.6–99.8%.','No pathogenic frameshift conclusion is supported for this allele.'),
'DST':('Lower priority; benign-leaning conflict','One uncertain and three likely-benign contributing submissions. MANE nomenclature is p.Glu4832Gly, whereas the original table selected p.Glu2209Gly on another transcript. Population frequency varies across datasets. A single heterozygous allele does not establish a recessive diagnosis.','Reconsider if phenotype fits the specific DST disorder and transcript; no second-allele investigation was completed here.'),
'ENPP1':('Deprioritize; common allele','ClinVar contains one likely-pathogenic submission versus eleven benign submissions. Listed population frequencies are approximately 19–36%, depending on dataset. This common allele is unsuitable as a rare, highly penetrant monogenic explanation on the evidence reviewed.','Do not infer diabetes, obesity, or a calcification disorder from the condition names attached to submissions.'),
'DLX6':('Unresolved; phenotype-dependent review','In-frame repeat insertion, not a frameshift. One uncertain and one likely-benign contributing submission. A submitter flags gnomAD frequency quality as unreliable; no dependable population frequency is established here. The repetitive sequence complicates alignment.','If clinically relevant, review in a genome browser and confirm the repeat allele before pursuing inheritance or functional evidence.'),
'STXBP1':('Lower priority; benign-leaning conflict','Synonymous p.Ile440= in MANE Select and Plus Clinical transcripts. One uncertain, one benign, and six likely-benign contributing submissions. ClinVar lists gnomAD frequencies around 0.006–0.008%. Synonymous status alone does not exclude splice effects, but no validated splice defect was established here.','Reconsider with compatible seizures/developmental history and variant-specific splice evidence; gene membership is insufficient.'),
'NOTCH1':('Unresolved; phenotype-dependent review','Missense p.Arg731Trp in MANE. Five uncertain, one benign, and one likely-benign contributing submissions. ClinVar lists low population frequencies (gnomAD approximately 0.004–0.016%, dataset-dependent); a submitter reports African/African-American frequency 0.038% in gnomAD v2.1.1. Predictions and interpretations disagree.','Clinical genetics review is most useful if there are relevant cardiac, vascular, scalp, or limb findings; this variant does not establish them.'),
'POLR3A':('Unresolved; phenotype-dependent review','A 5-prime UTR change (c.-35C>G), with one uncertain and one likely-benign contributing submission. ClinVar lists gnomAD frequencies approximately 0.12–0.14%. Only one heterozygous allele is present in this shortlist; a second disease-causing allele was not established.','Consider only in a relevant neurologic/dental/endocrine context and with inheritance-aware review; no recessive diagnosis follows from this row.'),
'PTEN':('Deprioritize; benign transcript-dependent flag','MANE places this deletion in the 5-prime UTR (c.-366del), despite a splice-site label on an alternative transcript. Four benign contributing submissions; the allele is extremely common in listed population resources.','This is not evidence of pathogenic PTEN loss of function or PTEN-related cancer/autism syndrome.'),
'SYCE1':('Deprioritize; likely-benign transcript-dependent flag','Likely benign, single submitter. The same position is upstream of the MANE transcript (c.-2798dup) rather than a MANE frameshift. A reliable population frequency is not supplied on the retrieved page.','The alternative-transcript frameshift label alone does not justify a disease association.'),
'CTR9':('Lower priority; benign-leaning conflict','Missense p.Gly1071Ser, with one uncertain and two likely-benign contributing submissions. ClinVar lists gnomAD frequencies roughly 0.046–0.105%, dataset-dependent. Prior reports in affected individuals do not establish causation; a submitter describes predictions as not damaging.','Reconsider only with relevant phenotype/inheritance evidence; no disease or tumor-risk conclusion is supported by this row.'),
'VDR':('Deprioritize; common allele','Intronic c.1025-49G>T, with one uncertain and three benign contributing submissions. Population frequency is approximately 51–56% across the listed resources.','No vitamin-D disorder or supplement recommendation follows from this variant.'),
'BRCA2':('Unresolved; incidental clinical review','Missense p.Asp651Tyr in MANE. Four uncertain and two likely-benign contributing submissions; none pathogenic in the contributing set. One submitter flags population-frequency quality as unreliable, so a rounded zero must not be interpreted as proven absence. The change lies in a region described by a submitter as a missense coldspot.','Discuss with clinical genetics if personal/family cancer history warrants it. This is not an actionable pathogenic BRCA2 finding; do not change cancer management on this VUS/conflict alone.'),
'MYO1E':('Deprioritize; common allele','Intronic/splice-region c.332+8G>A; one uncertain, four benign, and one likely-benign contributing submissions. The ALT allele, not the reference allele, occurs at approximately 97–99.8% in listed resources. Homozygosity therefore does not make it suspicious.','No monogenic kidney-disease inference is justified from this row.'),
'ZFHX3':('Lower priority; benign-leaning conflict','Synonymous p.Asn82= in MANE, intronic on another transcript. One uncertain and one benign contributing submission. Listed gnomAD frequencies are about 0.21–0.34%. No validated damaging or splice effect was established here.','Low priority absent additional variant-specific evidence; synonymous status alone is not a formal benign classification.'),
'POLR2A':('Deprioritize; benign reference/transcript discrepancy','Two benign contributing submissions. ClinVar maps this GRCh38 repeat-length difference to MANE c.5511_5512=, with no sequence alteration relative to that transcript. The corresponding MANE transcript was not present among the full-database ANN entries, so the local worst-effect frameshift must not be treated as established MANE loss of function.','Retain the reference/transcript discrepancy in the audit trail; the homozygous VCF call does not establish pathogenic POLR2A disruption.'),
'SETBP1':('Deprioritize; benign transcript-dependent flag','Five benign contributing germline submissions. The insertion is intronic in MANE but frameshifting in another transcript. Listed population frequencies are approximately 49–54%. Local ALT-supporting reads also carry a nearby flank substitution.','Do not conflate this common intronic MANE allele with pathogenic SETBP1 hotspot or loss-of-function variants.'),
'KATNAL2':('Deprioritize; benign, technically weak repeat call','ClinVar is benign from one submitter without assertion criteria, so review strength is limited. The insertion is intronic in MANE and frameshifting in some alternative transcripts. It occurs in a long A repeat; caller support is 5/23 ALT reads. The stricter local full-window check retains only one ALT-supporting read, on one strand.','Treat the exact repeat allele as technically uncertain pending visual/orthogonal confirmation if ever clinically relevant; no mosaicism inference is supported.')}
details=json.loads((OUT/'details.json').read_text());reads={r['gene']:r for r in json.loads((OUT/'read_support.json').read_text())}
mane=defaultdict(set)
with (ROOT/'data/reference/mane.features.tsv').open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['type']=='Transcript':mane[r['geneName']].add(r['id'].split('.')[0])
rows=[];paragraphs=[]
for d in details:
 r=d['candidate'];g=r['gene'];vid=IDS[g];s=(OUT/f'clinvar_{vid}.txt').read_text()
 a=s.find('Germline\nClassification');b=s.find('Somatic',a);summary=s[a:b]
 m=re.search(r'(?:Uncertain significance|Likely pathogenic|Benign|Likely benign) \(\d+\)(?:; [^\n]+)?',summary)
 assert m,g
 counts=m.group(0);category,note,nextstep=NOTES[g];url=f'https://www.ncbi.nlm.nih.gov/clinvar/variation/{vid}/'
 annotations=[a for a in d['annotations'] if a[6].split('.')[0] in mane[g]]
 manetext='; '.join(f'{a[6]}: {a[1]} {a[9]} {a[10]}'.strip() for a in annotations) or 'No matching MANE transcript ANN; see current ClinVar MANE mapping.'
 rr=reads[g];refn=rr.get('ref',0);altn=rr.get('alt',0)
 row={**r,'review_disposition':category,'clinvar_variation_id':vid,'current_contributing_classifications':counts,'mane_annotation':manetext,'review_reason':note,'suggested_followup':nextstep,'local_window_ref_reads':refn,'local_window_alt_reads':altn,'local_window_other_reads':rr.get('other',0),'local_window_alt_forward':rr.get('alt_forward',0),'local_window_alt_reverse':rr.get('alt_reverse',0),'public_source':url}
 rows.append(row)
 paragraphs.append(f'''### {g}: chr{r['chromosome']}:{r['position']} {r['ref']}>{r['alt']} (GRCh38, VCF representation)

**{category}.** {note}

Original caller: GT {r['genotype']}, DP {r['depth']}, GQ {r['genotype_quality']}, AD {r['allele_depths']}. MANE comparison: {manetext}

Local window support: REF {refn}, ALT {altn}, other {rr.get('other',0)}; ALT forward/reverse {rr.get('alt_forward',0)}/{rr.get('alt_reverse',0)}. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: {nextstep}

Source: [ClinVar variation {vid}]({url}); contributing submissions: {counts}.
''')
with (OUT/'571566.reviewed_candidates.tsv').open('w') as f:
 w=csv.DictWriter(f,rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
assert len(rows)==18 and len({r['gene'] for r in rows})==18
header='''# Review of the 18 shortlisted variants: sample 571566

Reviewed on 2026-10-20. This is a research triage of the original 18-row shortlist, not a clinical diagnosis or a complete exome interpretation. No phenotype, parental genotypes, segregation data, or detailed family history were supplied for this review.

None of these 18 variants is established here as a pathogenic explanation for autism or another disorder. Six entries have benign/likely-benign classifications and transcript/reference-context issues; three further entries are common alleles; five have benign-leaning or low-suspicion conflicting evidence; four remain unresolved for phenotype-dependent review. These are practical review groups, not new ACMG classifications.

The unresolved entries are DLX6, NOTCH1, POLR3A, and BRCA2. Their inclusion does not mean they are likely disease-causing. BRCA2 is an incidental cancer-gene uncertainty, not an established autism explanation. No treatment or screening change is justified by this review alone. A VUS should not guide clinical decisions by itself ([ACMG/AMP guidance](https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/)).

The original screen selected the worst predicted effect in any transcript. It also matched the word “pathogenic” inside “Conflicting classifications of pathogenicity”; that label often represents a VUS-versus-benign disagreement. This explains why the 18-row shortlist contained many weak leads. The original report and shortlist are retained; this review supplements them.

## Methods and limits

- Matched all 18 records by GRCh38 position, REF and ALT in the local normalized ClinVar VCF. ClinVar Variation IDs are distinct from Allele IDs; source links use the matched Variation IDs.
- Downloaded all 18 public ClinVar variant pages and cached their HTML/text locally. Submission counts below describe contributing submissions, not independent votes or a new formal classification. Conditions and submissions can differ, and a more recent submission is not automatically better evidence.
- Compared all full Ensembl 115 ANN entries with locally available MANE transcript IDs, ignoring transcript-version suffixes for matching. ClinVar MANE mappings provide an additional check where the full database lacks a transcript. Alternative-transcript effects are not automatically irrelevant, but need disease-relevant transcript evidence before a loss-of-function inference.
- Population evidence comes from the public ClinVar frequency tables and submitter narratives, plus the original limited 1000 Genomes annotation. ClinVar displays multiple dataset releases; ranges are descriptive, not ancestry-specific maximum frequencies. A separate comprehensive gnomAD v4.1 query, coverage/quality assessment, homozygote count review, and formal FAF/ACMG frequency analysis were not performed. Missing or rounded-zero frequency is not proof of rarity.
- Re-read the recalibrated BAM locally with samtools. Excluded unmapped, secondary, supplementary, QC-failed and duplicate reads; required mapping quality >=20, a read spanning the variant plus 25 bases on each side, and base quality >=20 throughout that window. Reconstructed CIGAR haplotypes and compared REF/ALT, allowing at most two flank-only substitutions while preserving the exact focal allele. Counts are intentionally more restrictive than GATK AD and can differ due to local realignment, adjacent variants, clipping and base qualities. Retained strand counts and other haplotypes in read_support.json.
- This was an automated local read-support check, not visual IGV inspection, independent validation, phase determination or mosaicism detection. In particular, the KATNAL2 homopolymer allele remains technically weak.
- This review does not search for a second recessive allele, assess the rest of the exome, resolve X/Y ploidy, or evaluate CNVs, structural variants, repeat expansions, mitochondrial disease or poorly covered regions.
- No sample BAM, FASTQ, or VCF was uploaded. Public evidence requests used known variant IDs and gene/variant search terms.

## All 18 candidates

'''
(OUT/'REVIEW.md').write_text(header+'\n'.join(paragraphs)+'''
## What to do next

1. Supply age, main developmental/neurologic symptoms, seizure history, hearing/vision findings, congenital findings and relevant family history. Review uncertain variants against those findings with a clinical geneticist/genetic counselor.
2. If a specific variant becomes clinically plausible, visually inspect its reads and confirm it with a clinically validated method before acting. Parental testing can help only when interpreted with the gene's inheritance mechanism and phenotype; do not order testing for every row indiscriminately.
3. If the question remains unexplained, a clinical review may need to extend beyond this SFARI-based shortlist, including appropriate CNV/repeat/sex-chromosome assessment and coverage review. This 18-variant review cannot provide a negative genetic diagnosis.

## Audit files

- 571566.reviewed_candidates.tsv: original fields plus MANE comparison, public classifications, technical support and review disposition.
- details.json: complete local annotations and caller INFO for these variants.
- read_support.json: conservative per-variant local haplotype/strand support.
- clinvar_<VariationID>.html/.txt: retrieved public records.
- provenance.json: input/public-page hashes and retrieval timestamps.
''')
files=[ROOT/'results/wgs/571566.sfari_priority_review.tsv',ROOT/'data/reference/mane.features.tsv',OUT/'details.json',OUT/'read_support.json',*sorted(OUT.glob('clinvar_*.html'))]
provenance={'built_utc':datetime.now(timezone.utc).isoformat(),'files':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'mtime_utc':datetime.fromtimestamp(p.stat().st_mtime,timezone.utc).isoformat()} for p in files]}
(OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
print('Saved review of',len(rows),'candidates to',OUT/'REVIEW.md')
