# Variant quality screen — sample 571566

Research quality selection, not a clinical validation or pathogenicity classification. Original VCFs are unchanged. Rejected records are retained for review, not declared false positives.

| Result | Records |
|---|---:|
| total | 97,174 |
| NEEDS_REVIEW | 40,645 |
| original_PASS | 73,378 |
| original_PASS_now_review | 16,849 |
| ACCEPTED | 56,529 |

## Rules

| Metric | Required |
|---|---|
| Original FILTER | PASS |
| QUAL / sample DP / sample GQ | ≥30 / ≥20 / ≥30 |
| QD / MQ | ≥2 / ≥40 (MQ extended to indels as a conservative local rule) |
| FS | ≤60 SNP; ≤200 other small alleles |
| SOR | ≤3 SNP; no extra indel SOR cutoff |
| ReadPosRankSum | ≥−8 SNP; ≥−20 indel, when available |
| MQRankSum | ≥−12.5 SNP, when available |
| Allelic depths | Two nonnegative AD entries, at least 5 ALT reads |
| Heterozygous ALT fraction | 0.25–0.75 |
| Homozygous/haploid ALT fraction | ≥0.90 |
| Existing independent SNV recount | At least 5 ALT reads at BQ20/MQ40 if a recount is available |
| Missing values | Required metrics fail; rank sums may be unavailable and are flagged |

These are conservative research thresholds. DP, GQ, allele balance and the indel MQ rule are local choices, not a clinically calibrated model. Low-depth genuine variants, mosaics and imbalanced calls can be excluded. No maximum-depth filter was imposed on exome data. Population frequency, predicted impact and ClinVar significance do not affect quality acceptance. Sex-chromosome ploidy is unresolved and flagged. Multi-ALT, symbolic and spanning-deletion placeholder records are retained for review. The existing independent recount covers only a small subset; accepted variants have not all had read-level alignment validation.

Site-filter basis: [GATK hard-filtering guidance](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants).

## Failure reasons (overlap possible)

| Reason | Records |
|---|---:|
| DP_below_20 | 35,750 |
| original_filter:LowGenotypeSupport | 15,909 |
| ALT_reads_below_5 | 14,602 |
| GQ_below_30 | 12,728 |
| heterozygous_balance_outside_0.25_0.75 | 5,438 |
| original_filter:SNPQuality | 4,619 |
| SOR_above_3 | 4,383 |
| MQ_below_40 | 3,149 |
| original_filter:SNPQuality;LowGenotypeSupport | 2,919 |
| QD_below_2 | 745 |
| original_filter:IndelQuality | 349 |
| homozygous_or_haploid_ALT_fraction_below_0.90 | 125 |
| unsupported_or_placeholder_allele | 110 |
| independent_ALT_MQ40_reads_below_5 | 2 |
| missing_or_invalid_AD | 2 |
| FS_above_60 | 1 |

## Findings from report genes

| Gene | GRCh38 allele | GT | DP | GQ | ALT fraction | Result | Reasons |
|---|---|---|---:|---:|---:|---|---|
| DNAH14 | 1:225240648 C>T | 0/1 | 34 | 99 | 0.44118 | ACCEPTED | — |
| WNT10A | 2:218889989 C>T | 1/1 | 66 | 99 | 0.98485 | ACCEPTED | — |
| DOLK;ENSG00000251184;NUP188;PHYHD1 | 9:128945903 AG>A | 0/1 | 90 | 99 | 0.47778 | ACCEPTED | — |
| VWF | 12:6034812 C>T | 0/1 | 87 | 99 | 0.44828 | ACCEPTED | — |
| ENSG00000258892;SLC38A6;TRMT5 | 14:60979582 CTATT>C | 0/1 | 78 | 99 | 0.53846 | ACCEPTED | — |
| FDXR;GRIN2C | 17:74863914 G>A | 1/1 | 36 | 67 | 0.97222 | ACCEPTED | — |
| ABCA7;CNN2;ENSG00000279753 | 19:1037646 AC>A | 0/1 | 75 | 99 | 0.24000 | NEEDS_REVIEW | heterozygous_balance_outside_0.25_0.75 |
| ABCA7;CNN2;ENSG00000279753 | 19:1037757 G>A | 0/1 | 70 | 99 | 0.21429 | NEEDS_REVIEW | heterozygous_balance_outside_0.25_0.75;independent_ALT_MQ40_reads_below_5 |
| ABCA7;CNN2;ENSG00000279753 | 19:1037767 G>A | 0/1 | 69 | 99 | 0.18841 | NEEDS_REVIEW | heterozygous_balance_outside_0.25_0.75;independent_ALT_MQ40_reads_below_5 |

Full per-record metrics, decisions and warnings: `variant_quality.tsv`. VCF subset FILTER columns retain their original meaning; the new screen decision is recorded in the TSV and by subset membership. The earlier clinical-style report remains a historical research draft; consult this quality addendum alongside it.
