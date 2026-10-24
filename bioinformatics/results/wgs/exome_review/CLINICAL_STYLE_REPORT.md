# Exome sequencing — clinical review draft

**Research findings requiring clinical confirmation. This document is not a validated diagnostic laboratory report.**

| Report information | Details |
|---|---|
| Sample identifier | 571566 |
| Report date | 20 October 2026 |
| Indication | Autism-related genetic investigation requested; clinical phenotype not provided |
| Analysis | Single-end whole-exome sequencing; small-variant review |
| Reference genome | GRCh38 |
| Patient name, date of birth, sex | Not provided |
| Family history / parental samples | Not provided / not evaluated |
| Clinical laboratory authorization | None; research review draft |

## Interpretation

Two apparently homozygous variants in **FDXR** and **WNT10A** are the principal findings for clinical confirmation. A heterozygous **VWF** variant is an additional finding. Public database classifications below describe the variants, not a diagnosis in this person. Clinical relevance cannot be established without confirmation and correlation with the person's features. **No molecular explanation for autism has been established by this review.**

## Principal findings

| Priority / gene | Transcript and variant | GRCh38 genomic allele | Observed genotype | Reported ClinVar classification | Interpretation for this sample |
|---|---|---|---|---|---|
| 1 — FDXR | NM_024417.5:c.1156C>T; p.Arg386Trp | chr17:74863914 G>A | Apparently homozygous | Pathogenic/likely pathogenic; multiple submitters, no conflicts; [520994](https://www.ncbi.nlm.nih.gov/clinvar/variation/520994/) | Potentially relevant to recessive FDXR-related disease with neurologic, optic or auditory manifestations. Confirm allele and genotype; phenotype correlation required. |
| 2 — WNT10A | NM_025216.3:c.382C>T; p.Arg128Ter | chr2:218889989 C>T | Apparently homozygous | Pathogenic/likely pathogenic in the retained review; [464181](https://www.ncbi.nlm.nih.gov/clinvar/variation/464181/) | Potentially relevant to WNT10A-related dental/ectodermal disease. Confirm allele and genotype; phenotype correlation required. |
| 3 — VWF | NM_000552.5:c.2561G>A; p.Arg854Gln | chr12:6034812 C>T | Heterozygous | Expert-panel pathogenic for recessive type 2N von Willebrand disease in the retained review; [296](https://www.ncbi.nlm.nih.gov/clinvar/variation/296/) | A single heterozygous allele does not establish type 2N disease. No second pathogenic coding allele was established; interpretation requires clinical review. |

Apparently homozygous calls do not establish two intact inherited copies. Copy number, allele dropout and segregation have not been assessed. No independent ACMG/AMP classification was performed. FDXR also appears as p.Arg392Trp on a different transcript; the transcript above specifies the reported nomenclature.

## Supporting sequence evidence

| Gene | Caller genotype quality | Caller REF / ALT reads | Independent REF / ALT reads | ALT forward / reverse | ALT mapping quality |
|---|---|---|---|---|---|
| FDXR | 67 | 1 / 35 | 1 / 34 | 18 / 16 | All 60 |
| WNT10A | 99 | 1 / 65 | 1 / 65 | 38 / 27 | All 60 |
| VWF | 99 | 48 / 39 | 48 / 39 | 18 / 21 | All 60 |

Independent counts are a computational recheck of the same BAM, requiring base and mapping quality at least 20, excluding duplicate and nonprimary/QC-failed reads. They are not independent laboratory confirmation. Differences from caller counts reflect differing filters and local realignment.

## Additional research findings

| Gene / variant | Observed genotype | Evidence status | Reporting interpretation |
|---|---|---|---|
| DNAH14 c.6574C>T; p.Gln2192Ter | Heterozygous | Single research pathogenic assertion; gene–disease validity and inheritance need review | Research lead; no second established pathogenic allele |
| TRMT5 c.312_315del; p.Ile105fs | Heterozygous | Conflicting pathogenic/likely pathogenic versus uncertain assertions | Unresolved; no established diagnosis or confirmed carrier interpretation |
| DOLK c.1400delC; p.Pro467fs | Heterozygous | Predicted frameshift; clinical classification and reliable frequency not established | Unclassified research lead |
| CNN2 p.Gly263Ser / p.Arg266Gln | Unresolved | Substantial alignment concerns; no ALT reads retained at mapping quality 40 | Excluded from principal findings |
| Other MANE HIGH candidates | Various | 128 candidates remain an unclassified research queue | Manual clinical review is incomplete; these are not reportable pathogenic findings |

The abbreviated HGVS descriptions in this additional-findings table are retained research annotations, not independently validated clinical nomenclature.

## Coverage and analysis scope

| Measure | Result |
|---|---|
| Annotated VCF records screened | 97,174 |
| Records passing original filters | 73,378 |
| Genes assessed using primary-chromosome MANE coding sequence plus 2 bp flanks | 19,229 |
| Genes with at least 90% of assessed bases at 20× | 17,055 |
| Genes below that threshold | 2,174, including 1,004 with no assessed bases at 20× |
| FDXR / WNT10A assessed bases at 20× | 100% / 100% |
| VWF assessed bases at 20× | 96.36% |
| SYNGAP1 assessed bases at 20× | 0.17%; major detection limitation |
| SHANK3 assessed bases at 20× | 87.62% in an earlier Ensembl supplement; absent from the primary MANE set |

Coverage uses mapping quality at least 20 with duplicate, secondary, supplementary, QC-failed and unmapped reads excluded, without a base-quality threshold. These are annotation-based measurements, not proof of capture-kit completeness. VCF record counts differ from the application's allele counts because spanning-deletion placeholder alleles are excluded there.

## Limitations

| Area | Limitation and implication |
|---|---|
| Clinical information | No phenotype, examination, pedigree or independent identity verification was available. BAM and VCF labels agree with 571566 only. |
| Confirmation | No independent sequencing, visual IGV review, parental segregation or phase determination was performed. |
| Variant classes | Copy-number changes, structural variants, repeat expansions and mitochondrial variants were not evaluated by this review. |
| Sensitivity | Poorly covered regions and variants outside the reviewed small-variant scope may be missed. A negative finding does not exclude a genetic cause. |
| Interpretation completeness | Full modern population annotation, ancestry-aware evaluation and manual review of all candidates are incomplete. The original 18-variant SFARI review was a restricted shortlist. |
| Clinical use | Findings require a clinically validated workflow and specialist interpretation before treatment or screening decisions. |

## Follow-up for the receiving clinician

| Step | Purpose |
|---|---|
| Confirm sample identity, exact alleles and genotypes using a clinically validated method | Establish whether the research findings are present and correctly genotyped |
| Review phenotype and family history with clinical genetics | Determine relevance of FDXR, WNT10A and VWF to this person |
| Select segregation and copy-number evaluation as clinically appropriate | Clarify inheritance and apparently homozygous calls |
| Assess remaining diagnostic needs, including poorly covered genes | Address the incomplete autism/developmental investigation |

## Evidence provenance

Underlying analysis and measurements: [Expanded review](REVIEW.md), [clinical handoff](CLINICAL_HANDOFF.md), original annotated VCF `../571566.annotated.vcf.gz`, retained read-support JSON and coverage tables. No new sequencing or reanalysis was performed to format this report. ClinVar classifications for WNT10A and VWF are carried forward from the retained review; FDXR was also checked against the public record during report preparation. This document is unsigned and has not been approved by a clinical laboratory.
