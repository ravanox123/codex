# Expanded exome review: 571566

Research review, 2026-10-20. This extends the original SFARI-only shortlist. It does not replace clinical confirmation or phenotype-based interpretation.

## Findings requiring clinical confirmation

| Priority | GRCh38 allele | MANE effect | Caller genotype / AD | Independent local read support (REF/ALT, BQ>=20, MQ>=20) | Evidence |
|---|---|---|---|---|---|
| 1 | chr17:74863914 G>A | FDXR c.1156C>T, p.Arg386Trp | 1/1; 1,35; GQ 67 | 1/34; ALT forward/reverse 18/16; all ALT MQ60 | [ClinVar 520994](https://www.ncbi.nlm.nih.gov/clinvar/variation/520994/) pathogenic/likely pathogenic, multiple submitters. |
| 2 | chr2:218889989 C>T | WNT10A c.382C>T, p.Arg128Ter | 1/1; 1,65; GQ 99 | 1/65; ALT forward/reverse 38/27; all ALT MQ60 | [ClinVar 464181](https://www.ncbi.nlm.nih.gov/clinvar/variation/464181/) pathogenic/likely pathogenic, multiple submitters. |
| 3 | chr12:6034812 C>T | VWF c.2561G>A, p.Arg854Gln | 0/1; 48,39; GQ 99 | 48/39; ALT forward/reverse 18/21; all ALT MQ60 | [ClinVar 296](https://www.ncbi.nlm.nih.gov/clinvar/variation/296/) expert-panel pathogenic for recessive type 2N von Willebrand disease. A single heterozygous allele does not establish that diagnosis. |

FDXR is associated with recessive mitochondrial disease with neurologic, optic and auditory manifestations. WNT10A is associated with dental/ectodermal disorders. These two apparently homozygous calls merit clinical genetics review and confirmation of the exact alleles, sample identity and genotype. These findings do not establish autism or show that either explains the donor's symptoms. No phenotype or family data were available. [FDXR evidence](https://www.ncbi.nlm.nih.gov/clinvar/variation/520994/), [WNT10A evidence](https://www.ncbi.nlm.nih.gov/clinvar/variation/464181/).

## Other reviewed leads

- DNAH14 c.6574C>T, p.Gln2192Ter: heterozygous; caller AD 19,15 and independent support 20 REF/15 ALT. One research submitter classifies it pathogenic. Its gene/disease relationship, assertion and inheritance need specialist review; no second pathogenic allele was established. [ClinVar 4880149](https://www.ncbi.nlm.nih.gov/clinvar/variation/4880149/).
- TRMT5 c.312_315del, p.Ile105fs: heterozygous; caller AD 36,42; conservative full-window support 21 REF/24 ALT. ClinVar has pathogenic/likely-pathogenic versus uncertain disagreement. No second qualifying coding allele was identified in the focused table; this is not a complete exclusion of another allele. [ClinVar 372246](https://www.ncbi.nlm.nih.gov/clinvar/variation/372246/).
- DOLK c.1400delC, p.Pro467fs: heterozygous predicted frameshift; conservative local support 23 REF/26 ALT. No clinical classification or dependable population frequency was established. Do not treat a computational frameshift alone as a pathogenic diagnosis or confirmed carrier result.
- CNN2 p.Gly263Ser and p.Arg266Gln: both have single-submitter, no-criteria pathogenic assertions. The ALT-supporting reads have mean MQ about 26 and median 13 edit differences versus median zero for REF reads; neither ALT survives MQ>=40. These calls have substantial alignment concerns and are excluded from the main findings. [ClinVar 1696852](https://www.ncbi.nlm.nih.gov/clinvar/variation/1696852/), [ClinVar 1696851](https://www.ncbi.nlm.nih.gov/clinvar/variation/1696851/).
- CDKN2B rs1063192: common allele; the retrieved current page reports a protective germline assertion rather than the historical local likely-pathogenic/protective aggregate. Not a rare monogenic finding. [ClinVar 812639](https://www.ncbi.nlm.nih.gov/clinvar/variation/812639/).
- PRKAA2 stop, CNTN5 stop, SMARCA5 indel, SUZ12 splice-region indel and RBMX indel remain computational/technical leads, not pathogenic classifications. Several have weak allele balance, repeat/nearby-indel issues, missing population evidence or unresolved sex-chromosome ploidy. The conservative window check retained no qualifying SMARCA5 reads. Population/gene-disease validation was not completed for these leads.

## Coverage

One sequential samtools bedcov pass measured primary-chromosome MANE CDS plus 2bp flanks in 19,229 genes. Mapping quality >=20; duplicate, secondary, supplementary, QC-failed and unmapped reads excluded; no base-quality filter. 17,055 genes had >=90% of these bases at 20x; 2,174 were below that threshold, including 1,004 with no bases at 20x. These measurements are annotation-based, not capture-kit completeness or CNV calls.

FDXR and WNT10A measured coding coverage reached 100% at 20x; VWF 96.36%. Important gaps remain: SYNGAP1 measured only 0.17% at 20x. SHANK3 is absent from the primary MANE set; the earlier Ensembl supplement measured 87.62% at 20x. Lack of calls cannot exclude variants in poorly covered regions.

## Scope and audit

All 97,174 annotated variants were screened; 73,378 passed the original filters. Among passing records, 545 have a HIGH prediction on some transcript and 228 on a matched MANE transcript. Seven records had explicit aggregate pathogenic/likely-pathogenic tokens in the local ClinVar annotation. All seven were checked against current public records and independently re-counted in the BAM. This is broader than the original 18-row SFARI review.

The 128 other MANE HIGH candidates without known common/benign labels remain an unclassified research queue. They have not all undergone manual clinical review, ancestry-aware frequency checks or gene-disease validation. All candidates, including filtered and alternative-transcript calls, are retained in `all_HIGH_or_pathogenic_assertion.tsv`; `mane_HIGH_candidates.tsv` preserves the queue. Nearby indels are flagged because joint haplotypes can restore reading frame; proximity alone does not establish phase. The 210-row `focused_gene_variants.tsv` retains passing and failed calls in 13 selected genes for further review.

Read checks used local CIGAR reconstruction and focal SNV counting on the recalibrated BAM, excluding duplicates and nonprimary/QC-failed reads. Counts may differ from GATK AD because filters and local realignment differ. No independent sequencing, visual IGV review, copy-number analysis, repeat-expansion testing, phase determination, full modern population dataset annotation, or parental segregation testing was performed. Sample labels in BAM and VCF match 571566, but that is not independent identity verification.

Public ClinVar HTML/text snapshots are cached in `../candidate_review/clinvar_<VariationID>.*`. No BAM, FASTQ or VCF was uploaded. The original [18-candidate review](../candidate_review/REVIEW.md) remains available; its conclusions apply only to those 18 rows.

## Clinical next step

Take the FDXR and WNT10A findings, VWF result and this report to clinical genetics for confirmation and phenotype review. Inheritance testing and clinical assessment can then be selected appropriately. These research files should not independently drive treatment or screening changes.
