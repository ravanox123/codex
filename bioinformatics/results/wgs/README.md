# 571566: single-end WES variant analysis

Analysis was started on 2026-10-16. This directory contains research outputs;
it does not establish or exclude an autism diagnosis.
The user confirmed this sample is WES. The existing directory name `wgs` is
retained so that running jobs and their output paths remain valid.

## Current progress

Read `STATUS.txt` for the pipeline's current stage and `logs/pipeline.log` for
stage history. Detailed logs are in `logs/`. `REPORT.md` and `summary.json`
are generated only after calling, annotation and screening finish successfully.

The analysis runs through `scripts/run_wgs.sh`, sequentially with a 2 GB Java
heap (3 GB for full-transcript annotation) and two HaplotypeCaller native threads.
Do not start another copy while the existing run is active. If a stage fails,
its log identifies the error and the status file records failure. After fixing
the issue, rerun the script; completed stages are skipped using completion files.
Completion files refer to this exact configuration and must not be reused for a
different sample or changed calling/filtering settings.

## Completed initial checks

- Original BAM passes samtools quickcheck and its index is readable.
- 46,979,864 primary reads, all single-end, 150 bases each.
- 46,818,523 primary reads mapped (99.66%).
- 12,582,587 mapped primary reads marked as duplicates (26.8752%).
- Median per-gene mean coverage across the measured SFARI coding intervals is
  62.54x after duplicate and mapping-quality filtering; 1,174 of 1,276 measured
  genes have at least 90% of their coding/flank bases at 20x.
- The overall aligned-base count is only about 2x averaged over the reference.
  Coding coverage is strongly enriched, as expected for WES. Whole-genome average
  depth does not measure exome completeness.
- After duplicate, mapping-quality and base-quality filtering, nuclear primary
  mean depth is 1.38x; 9.59% of bases reach 1x, 2.99% reach 10x and 2.28% reach 20x.
- The coverage table is `qc/571566.sfari_gene_coverage.tsv`; missing MANE coding
  genes are listed separately. This table uses mapping quality >=20, excludes
  duplicate/secondary/supplementary reads, and has no base-quality cutoff.
- Unrecalibrated small-region caller/filter/annotation files prefixed `qc/pilot`
  are technical smoke tests, not the final callset and not an interpretation.

## Scope and planned outputs

The original BAM is retained. A local copy of the exact GRCh38 FASTA is used.
Recalibration resources and ClinVar were matched to NCBI accession chromosome
names, checked against reference alleles, normalized, sorted and indexed.
The input lacks platform metadata; GATK's explicit UNKNOWN fallback preserves
that uncertainty during recalibration instead of assigning an instrument.

The script calls SNPs and small indels within merged Ensembl 115 exon intervals
plus 100bp flanks on the nuclear primary chromosomes. The exact capture-kit BED
was unavailable, so these are approximate annotation-derived WES targets;
kit-specific on-target rates and completeness cannot yet be measured.
It then
labels low-quality calls, annotates ClinVar and limited 1000 Genomes frequencies,
then predicts Ensembl 115 transcript effects with SnpEff and screens SFARI-listed genes.
Initial MANE coding coverage is supplemented with full Ensembl transcript unions
for genes missing on primary chromosomes in the MANE set, including checking SHANK3.
X/Y ploidy is unconfirmed: diploid X/Y calls are provisional and excluded from
the priority review table. Mitochondrial variants, alternate-contig independent
calls, structural variants, CNVs and repeat expansions are outside this workflow.

Expected final files:

- `571566.annotated.vcf.gz`: annotated calls with quality labels.
- `571566.pass.vcf.gz`: passing calls with DP>=10 and GQ>=20.
- `571566.sfari_all.tsv`: SFARI gene-overlapping variants and support notes.
- `571566.sfari_priority_review.tsv`: candidates requiring manual evidence review.
- `REPORT.md`: method, coverage, counts and interpretation limits.

Absence from the limited frequency resource does not establish rarity. SFARI
membership does not establish variant pathogenicity. ClinVar assertions may
conflict and may concern another condition. Inheritance, phenotype, read-level
review and clinical confirmation remain necessary for meaningful findings.

## Public resources

Only public reference resources were downloaded; sample genomic data was not uploaded.

- NCBI GRCh38 assembly report:
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_report.txt
- Broad hg38 known SNPs and Mills indels:
  https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/
- ClinVar GRCh38:
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
- SFARI current gene download:
  https://gene.sfari.org/wp-content/themes/sfari-gene/utilities/download-csv.php?api-endpoint=genes
- SnpEff GRCh38 MANE v1.5 Ensembl database:
  https://snpeff-public.s3.amazonaws.com/databases/v5_4/snpEff_v5_4_GRCh38.mane.1.5.ensembl.zip
- SnpEff full GRCh38 Ensembl 115 database:
  https://snpeff-public.s3.amazonaws.com/databases/v5_4/snpEff_v5_4_GRCh38.115.zip

Installed tool versions: GATK 4.7.0.0, samtools 1.24, bcftools 1.24,
SnpEff 5.4c. Downloaded resource copies and chromosome mappings are retained in
`data/reference/` and `data/snpeff/`.

## Sequential annotation chunks

Annotation now runs through `scripts/annotate_chunks.sh`, also invoked by
`run_wgs.sh`. It processes 10 Mb coordinate windows sequentially with a 3 GB
Java heap, one JVM processor and serial garbage collection. Full Ensembl 115
transcript context remains available; chromosome coordinates are preserved.
SnpEff still loads its database for each nonempty window, so chunking adds time
and does not eliminate database memory requirements.

Restart with `bash scripts/annotate_chunks.sh`. Valid completed windows in
`annotation_chunks/` are skipped. Input fingerprints prevent reuse after input
changes, and a lock prevents simultaneous chunk runners. Variants belong to a
window by their starting position; final merging checks both record count and
ordered chromosome/position/REF/ALT against the filtered input. Progress is in
`STATUS.txt`, `logs/annotation_chunks.log`, and per-window SnpEff logs. The final
VCF and report are published only after all windows succeed.
