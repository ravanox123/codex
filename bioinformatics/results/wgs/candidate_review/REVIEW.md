# Review of the 18 shortlisted variants: sample 571566

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

### USH2A: chr1:215671031 C>T (GRCh38, VCF representation)

**Lower priority; benign-leaning conflict.** Most contributing submissions are benign/likely benign, with one uncertain. ClinVar lists gnomAD frequencies around 0.42–0.68%, varying by dataset. This heterozygous change does not establish a recessive USH2A diagnosis.

Original caller: GT 0/1, DP 53, GQ 99, AD 26,27. MANE comparison: ENST00000307340.8: missense_variant c.14074G>A p.Gly4692Arg

Local window support: REF 20, ALT 18, other 0; ALT forward/reverse 12/6. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: Reconsider only with relevant hearing/retinal findings and a full inheritance-aware analysis.

Source: [ClinVar variation 48427](https://www.ncbi.nlm.nih.gov/clinvar/variation/48427/); contributing submissions: Uncertain significance (1); Benign (7); Likely benign (6).

### EXOC6B: chr2:72465143 A>AAG (GRCh38, VCF representation)

**Deprioritize; benign transcript-dependent flag.** The MANE transcript places the duplication in an intron (c.1980+15_1980+16dup), rather than the alternative-transcript frameshift selected by the screen. ClinVar is benign; population datasets listed there place the allele near 99.6–99.8%.

Original caller: GT 1/1, DP 28, GQ 84, AD 0,28. MANE comparison: ENST00000272427.11: intron_variant c.1980+15_1980+16dupCT

Local window support: REF 0, ALT 15, other 0; ALT forward/reverse 12/3. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: No pathogenic frameshift conclusion is supported for this allele.

Source: [ClinVar variation 1601386](https://www.ncbi.nlm.nih.gov/clinvar/variation/1601386/); contributing submissions: Benign (1).

### DST: chr6:56557464 T>C (GRCh38, VCF representation)

**Lower priority; benign-leaning conflict.** One uncertain and three likely-benign contributing submissions. MANE nomenclature is p.Glu4832Gly, whereas the original table selected p.Glu2209Gly on another transcript. Population frequency varies across datasets. A single heterozygous allele does not establish a recessive diagnosis.

Original caller: GT 0/1, DP 51, GQ 99, AD 24,27. MANE comparison: ENST00000680361.1: missense_variant c.14495A>G p.Glu4832Gly

Local window support: REF 11, ALT 18, other 0; ALT forward/reverse 9/9. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: Reconsider if phenotype fits the specific DST disorder and transcript; no second-allele investigation was completed here.

Source: [ClinVar variation 809956](https://www.ncbi.nlm.nih.gov/clinvar/variation/809956/); contributing submissions: Uncertain significance (1); Likely benign (3).

### ENPP1: chr6:131851228 A>C (GRCh38, VCF representation)

**Deprioritize; common allele.** ClinVar contains one likely-pathogenic submission versus eleven benign submissions. Listed population frequencies are approximately 19–36%, depending on dataset. This common allele is unsuitable as a rare, highly penetrant monogenic explanation on the evidence reviewed.

Original caller: GT 0/1, DP 53, GQ 99, AD 24,29. MANE comparison: ENST00000647893.1: missense_variant c.517A>C p.Lys173Gln

Local window support: REF 18, ALT 14, other 0; ALT forward/reverse 8/6. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: Do not infer diabetes, obesity, or a calcification disorder from the condition names attached to submissions.

Source: [ClinVar variation 13589](https://www.ncbi.nlm.nih.gov/clinvar/variation/13589/); contributing submissions: Likely pathogenic (1); Benign (11).

### DLX6: chr7:97006108 A>AGCCGCC (GRCh38, VCF representation)

**Unresolved; phenotype-dependent review.** In-frame repeat insertion, not a frameshift. One uncertain and one likely-benign contributing submission. A submitter flags gnomAD frequency quality as unreliable; no dependable population frequency is established here. The repetitive sequence complicates alignment.

Original caller: GT 0/1, DP 28, GQ 99, AD 12,16. MANE comparison: ENST00000518156.3: disruptive_inframe_insertion c.155_160dupCGCCGC p.Pro52_Pro53dup

Local window support: REF 11, ALT 10, other 1; ALT forward/reverse 6/4. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: If clinically relevant, review in a genome browser and confirm the repeat allele before pursuing inheritance or functional evidence.

Source: [ClinVar variation 1385415](https://www.ncbi.nlm.nih.gov/clinvar/variation/1385415/); contributing submissions: Uncertain significance (1); Likely benign (1).

### STXBP1: chr9:127676714 C>T (GRCh38, VCF representation)

**Lower priority; benign-leaning conflict.** Synonymous p.Ile440= in MANE Select and Plus Clinical transcripts. One uncertain, one benign, and six likely-benign contributing submissions. ClinVar lists gnomAD frequencies around 0.006–0.008%. Synonymous status alone does not exclude splice effects, but no validated splice defect was established here.

Original caller: GT 0/1, DP 88, GQ 99, AD 39,49. MANE comparison: ENST00000373302.8: synonymous_variant c.1320C>T p.Ile440Ile; ENST00000373299.5: synonymous_variant c.1320C>T p.Ile440Ile

Local window support: REF 29, ALT 34, other 0; ALT forward/reverse 14/20. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: Reconsider with compatible seizures/developmental history and variant-specific splice evidence; gene membership is insufficient.

Source: [ClinVar variation 365052](https://www.ncbi.nlm.nih.gov/clinvar/variation/365052/); contributing submissions: Uncertain significance (1); Benign (1); Likely benign (6).

### NOTCH1: chr9:136514526 G>A (GRCh38, VCF representation)

**Unresolved; phenotype-dependent review.** Missense p.Arg731Trp in MANE. Five uncertain, one benign, and one likely-benign contributing submissions. ClinVar lists low population frequencies (gnomAD approximately 0.004–0.016%, dataset-dependent); a submitter reports African/African-American frequency 0.038% in gnomAD v2.1.1. Predictions and interpretations disagree.

Original caller: GT 0/1, DP 77, GQ 99, AD 41,36. MANE comparison: ENST00000651671.1: missense_variant c.2191C>T p.Arg731Trp

Local window support: REF 29, ALT 27, other 0; ALT forward/reverse 16/11. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: Clinical genetics review is most useful if there are relevant cardiac, vascular, scalp, or limb findings; this variant does not establish them.

Source: [ClinVar variation 477890](https://www.ncbi.nlm.nih.gov/clinvar/variation/477890/); contributing submissions: Uncertain significance (5); Benign (1); Likely benign (1).

### POLR3A: chr10:78029442 G>C (GRCh38, VCF representation)

**Unresolved; phenotype-dependent review.** A 5-prime UTR change (c.-35C>G), with one uncertain and one likely-benign contributing submission. ClinVar lists gnomAD frequencies approximately 0.12–0.14%. Only one heterozygous allele is present in this shortlist; a second disease-causing allele was not established.

Original caller: GT 0/1, DP 26, GQ 99, AD 13,13. MANE comparison: ENST00000372371.8: 5_prime_UTR_variant c.-35C>G

Local window support: REF 8, ALT 7, other 0; ALT forward/reverse 1/6. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: Consider only in a relevant neurologic/dental/endocrine context and with inheritance-aware review; no recessive diagnosis follows from this row.

Source: [ClinVar variation 368925](https://www.ncbi.nlm.nih.gov/clinvar/variation/368925/); contributing submissions: Uncertain significance (1); Likely benign (1).

### PTEN: chr10:87864103 CT>C (GRCh38, VCF representation)

**Deprioritize; benign transcript-dependent flag.** MANE places this deletion in the 5-prime UTR (c.-366del), despite a splice-site label on an alternative transcript. Four benign contributing submissions; the allele is extremely common in listed population resources.

Original caller: GT 1/1, DP 26, GQ 78, AD 0,26. MANE comparison: ENST00000371953.8: 5_prime_UTR_variant c.-366delT

Local window support: REF 0, ALT 11, other 0; ALT forward/reverse 10/1. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: This is not evidence of pathogenic PTEN loss of function or PTEN-related cancer/autism syndrome.

Source: [ClinVar variation 440216](https://www.ncbi.nlm.nih.gov/clinvar/variation/440216/); contributing submissions: Benign (4).

### SYCE1: chr10:133568326 G>GC (GRCh38, VCF representation)

**Deprioritize; likely-benign transcript-dependent flag.** Likely benign, single submitter. The same position is upstream of the MANE transcript (c.-2798dup) rather than a MANE frameshift. A reliable population frequency is not supplied on the retrieved page.

Original caller: GT 0/1, DP 20, GQ 99, AD 9,11. MANE comparison: ENST00000343131.7: upstream_gene_variant c.-2798dupG

Local window support: REF 5, ALT 5, other 0; ALT forward/reverse 1/4. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: The alternative-transcript frameshift label alone does not justify a disease association.

Source: [ClinVar variation 2672423](https://www.ncbi.nlm.nih.gov/clinvar/variation/2672423/); contributing submissions: Likely benign (1).

### CTR9: chr11:10778794 G>A (GRCh38, VCF representation)

**Lower priority; benign-leaning conflict.** Missense p.Gly1071Ser, with one uncertain and two likely-benign contributing submissions. ClinVar lists gnomAD frequencies roughly 0.046–0.105%, dataset-dependent. Prior reports in affected individuals do not establish causation; a submitter describes predictions as not damaging.

Original caller: GT 0/1, DP 63, GQ 99, AD 26,37. MANE comparison: ENST00000361367.7: missense_variant c.3211G>A p.Gly1071Ser

Local window support: REF 17, ALT 23, other 0; ALT forward/reverse 12/11. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: Reconsider only with relevant phenotype/inheritance evidence; no disease or tumor-risk conclusion is supported by this row.

Source: [ClinVar variation 1091039](https://www.ncbi.nlm.nih.gov/clinvar/variation/1091039/); contributing submissions: Uncertain significance (1); Likely benign (2).

### VDR: chr12:47845054 C>A (GRCh38, VCF representation)

**Deprioritize; common allele.** Intronic c.1025-49G>T, with one uncertain and three benign contributing submissions. Population frequency is approximately 51–56% across the listed resources.

Original caller: GT 0/1, DP 49, GQ 99, AD 25,24. MANE comparison: ENST00000549336.6: intron_variant c.1025-49G>T

Local window support: REF 9, ALT 8, other 0; ALT forward/reverse 0/8. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: No vitamin-D disorder or supplement recommendation follows from this variant.

Source: [ClinVar variation 1252227](https://www.ncbi.nlm.nih.gov/clinvar/variation/1252227/); contributing submissions: Uncertain significance (1); Benign (3).

### BRCA2: chr13:32336306 G>T (GRCh38, VCF representation)

**Unresolved; incidental clinical review.** Missense p.Asp651Tyr in MANE. Four uncertain and two likely-benign contributing submissions; none pathogenic in the contributing set. One submitter flags population-frequency quality as unreliable, so a rounded zero must not be interpreted as proven absence. The change lies in a region described by a submitter as a missense coldspot.

Original caller: GT 0/1, DP 69, GQ 99, AD 39,30. MANE comparison: ENST00000380152.8: missense_variant c.1951G>T p.Asp651Tyr

Local window support: REF 23, ALT 21, other 0; ALT forward/reverse 16/5. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: Discuss with clinical genetics if personal/family cancer history warrants it. This is not an actionable pathogenic BRCA2 finding; do not change cancer management on this VUS/conflict alone.

Source: [ClinVar variation 51231](https://www.ncbi.nlm.nih.gov/clinvar/variation/51231/); contributing submissions: Uncertain significance (4); Likely benign (2).

### MYO1E: chr15:59256276 C>T (GRCh38, VCF representation)

**Deprioritize; common allele.** Intronic/splice-region c.332+8G>A; one uncertain, four benign, and one likely-benign contributing submissions. The ALT allele, not the reference allele, occurs at approximately 97–99.8% in listed resources. Homozygosity therefore does not make it suspicious.

Original caller: GT 1/1, DP 35, GQ 99, AD 0,35. MANE comparison: ENST00000288235.9: splice_region_variant&intron_variant c.332+8G>A

Local window support: REF 0, ALT 24, other 0; ALT forward/reverse 13/11. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: No monogenic kidney-disease inference is justified from this row.

Source: [ClinVar variation 983002](https://www.ncbi.nlm.nih.gov/clinvar/variation/983002/); contributing submissions: Uncertain significance (1); Benign (4); Likely benign (1).

### ZFHX3: chr16:72959900 G>A (GRCh38, VCF representation)

**Lower priority; benign-leaning conflict.** Synonymous p.Asn82= in MANE, intronic on another transcript. One uncertain and one benign contributing submission. Listed gnomAD frequencies are about 0.21–0.34%. No validated damaging or splice effect was established here.

Original caller: GT 0/1, DP 62, GQ 99, AD 40,22. MANE comparison: ENST00000268489.10: synonymous_variant c.246C>T p.Asn82Asn

Local window support: REF 20, ALT 14, other 0; ALT forward/reverse 11/3. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: Low priority absent additional variant-specific evidence; synonymous status alone is not a formal benign classification.

Source: [ClinVar variation 808107](https://www.ncbi.nlm.nih.gov/clinvar/variation/808107/); contributing submissions: Uncertain significance (1); Benign (1).

### POLR2A: chr17:7513774 CCA>C (GRCh38, VCF representation)

**Deprioritize; benign reference/transcript discrepancy.** Two benign contributing submissions. ClinVar maps this GRCh38 repeat-length difference to MANE c.5511_5512=, with no sequence alteration relative to that transcript. The corresponding MANE transcript was not present among the full-database ANN entries, so the local worst-effect frameshift must not be treated as established MANE loss of function.

Original caller: GT 1/1, DP 104, GQ 99, AD 0,104. MANE comparison: No matching MANE transcript ANN; see current ClinVar MANE mapping.

Local window support: REF 0, ALT 65, other 0; ALT forward/reverse 27/38. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: Retain the reference/transcript discrepancy in the audit trail; the homozygous VCF call does not establish pathogenic POLR2A disruption.

Source: [ClinVar variation 1686168](https://www.ncbi.nlm.nih.gov/clinvar/variation/1686168/); contributing submissions: Benign (2).

### SETBP1: chr18:44876705 C>CTCTT (GRCh38, VCF representation)

**Deprioritize; benign transcript-dependent flag.** Five benign contributing germline submissions. The insertion is intronic in MANE but frameshifting in another transcript. Listed population frequencies are approximately 49–54%. Local ALT-supporting reads also carry a nearby flank substitution.

Original caller: GT 0/1, DP 87, GQ 99, AD 47,40. MANE comparison: ENST00000649279.2: intron_variant c.540+7422_540+7423insTCTT

Local window support: REF 26, ALT 25, other 0; ALT forward/reverse 11/14. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: Do not conflate this common intronic MANE allele with pathogenic SETBP1 hotspot or loss-of-function variants.

Source: [ClinVar variation 326725](https://www.ncbi.nlm.nih.gov/clinvar/variation/326725/); contributing submissions: Benign (5).

### KATNAL2: chr18:47098221 C>CA (GRCh38, VCF representation)

**Deprioritize; benign, technically weak repeat call.** ClinVar is benign from one submitter without assertion criteria, so review strength is limited. The insertion is intronic in MANE and frameshifting in some alternative transcripts. It occurs in a long A repeat; caller support is 5/23 ALT reads. The stricter local full-window check retains only one ALT-supporting read, on one strand.

Original caller: GT 0/1, DP 23, GQ 63, AD 18,5. MANE comparison: ENST00000683218.1: intron_variant c.1212-1009dupA

Local window support: REF 7, ALT 1, other 1; ALT forward/reverse 1/0. This is a conservative technical cross-check, not orthogonal confirmation.

Follow-up: Treat the exact repeat allele as technically uncertain pending visual/orthogonal confirmation if ever clinically relevant; no mosaicism inference is supported.

Source: [ClinVar variation 3060277](https://www.ncbi.nlm.nih.gov/clinvar/variation/3060277/); contributing submissions: Benign (1).

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
