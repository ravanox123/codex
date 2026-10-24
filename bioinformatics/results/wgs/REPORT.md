# Single-end WES research analysis: 571566

This is a candidate-variant screen, not an autism diagnosis or a clinical genetic report.
Gene membership and computational effect predictions do not establish pathogenicity.
The sample is WES: low whole-genome coverage is expected and is not an exome-quality
failure. Coding-region coverage is the relevant measure. The exact capture-kit BED
was unavailable; annotation-derived exon targets cannot establish kit-specific
capture completeness or on-target enrichment.

## Coverage

- Nuclear primary bases at >= 1x: 9.59%.
- Mean nuclear primary depth after quality/duplicate filtering: 1.38x.
- Nuclear primary bases at >= 10x: 2.99%.
- Nuclear primary bases at >= 20x: 2.28%.
- SFARI genes with MANE coding coverage measurements: 1278.
- Genes with >=90% of measured coding/flank bases at >=20x: 1175.

Coverage excludes duplicate, secondary, supplementary, QC-failed and unmapped reads;
minimum mapping and base quality are both 20. Shallow or uncovered regions cannot
support negative findings. Single-end duplicate marking can also label independent
reads with identical alignment starts; optical duplicate counts are unavailable.

## Variant counts

- all_variants: 97,174
- pass_variants: 73,378
- sfari_gene_variants: 9,755
- sfari_supported_autosomal_variants: 3,706
- priority_review_variants: 18

## Files

- `571566.annotated.vcf.gz`: all annotated calls, including calls with failure labels.
- `571566.pass.vcf.gz`: calls passing the exploratory variant and genotype filters.
- `571566.sfari_all.tsv`: variants overlapping SFARI genes, with support and evidence notes.
- `571566.sfari_priority_review.tsv`: supported autosomal candidates with ClinVar
  pathogenic-related assertions (including conflicting assertions) or predicted HIGH
  effect without a known common frequency. These are review candidates, not diagnoses.
- `qc/571566.sfari_gene_coverage.tsv`: per-gene mean depth and fractions at 1x/10x/20x,
  across unioned MANE coding intervals with 2bp flanks, supplemented with Ensembl 115
  transcript unions for genes missing on primary chromosomes in the MANE set.
  Per-gene transcript sources are in `qc/sfari_gene_coverage_sources.tsv`.
  This separate bedcov measurement
  requires mapping quality >=20 and excludes duplicates/secondary/supplementary reads;
  it does not impose a base-quality threshold. Genes missing from this transcript set
  are listed in `qc/sfari_genes_without_coding_intervals.txt`.

## Method and limits

Original GRCh38 FASTA was copied locally and indexed. Known SNP/indel and ClinVar
resources were mapped using the NCBI assembly report, restricted to matching contigs,
reference-allele checked, normalized and sorted. Resource logs record exclusions.
GATK MarkDuplicates, BaseRecalibrator, ApplyBQSR and HaplotypeCaller were used.
The input lacks sequencing-platform metadata. Recalibration uses the explicit
UNKNOWN platform fallback (not an invented instrument label); its read-group,
quality, context and cycle covariates were smoke-tested. The actual sequencing
platform and the suitability of recalibration for it remain unconfirmed.
Calling is restricted to merged Ensembl 115 exon regions plus 100bp flanks on
the 24 nuclear primary chromosomes, as approximate WES targets. This includes
annotated UTR and noncoding exons; it is not the manufacturer's capture design.
Intergenic and deep-intronic variants are not comprehensively assessed.
Alternate contigs and mitochondrial
variants are not independently interpreted. X/Y calls use provisional diploid ploidy
and are excluded from the priority table pending appropriate ploidy information.

SNP failure thresholds: QUAL<30, QD<2, FS>60, SOR>3, MQ<40,
MQRankSum<-12.5 or ReadPosRankSum<-8. Indel failure thresholds: QUAL<30,
QD<2, FS>200 or ReadPosRankSum<-20. Missing DP/GQ, DP<10 or GQ<20
receive LowGenotypeSupport. Missing rank-sum annotations do not alone fail a site.
These are exploratory hard filters, not a clinically validated classifier.

SnpEff GRCh38 Ensembl 115 predicts transcript effects; a worst-effect annotation
in one transcript does not establish a clinically relevant loss of function.
ClinVar assertions can conflict or refer to conditions other than autism. The limited
1000 Genomes Phase 1 SNP resource supplies POP_AF when matched; missing frequency
does not establish rarity, and comprehensive modern population-frequency review is
still required. SFARI scores describe gene-level research evidence, not individual
variant pathogenicity. Inheritance, phenotype, segregation, read inspection and clinical
confirmation remain necessary for potentially significant findings.

Copy-number/structural variants, repeat expansions, polygenic scores, mosaicism and
comprehensive noncoding interpretation are not covered by this small-variant workflow.
Low coverage and single-end data further limit detection. No result rules out autism.
Public resources were downloaded; sample genomic data was not uploaded.

## Expanded exome review

A broader review beyond SFARI genes is available in [exome_review/REVIEW.md](exome_review/REVIEW.md), with [a clinical handoff](exome_review/CLINICAL_HANDOFF.md). It identifies additional FDXR and WNT10A calls requiring clinical confirmation; the original 18-row shortlist was not an exhaustive exome review.
