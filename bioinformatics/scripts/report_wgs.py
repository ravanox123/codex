"""Summarize a local research callset without assigning clinical diagnoses."""
import csv
import gzip
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/wgs'
REF = ROOT / 'data/reference'
genes = {}
ensembl = {}
with (REF / 'sfari-genes.csv').open() as f:
    for row in csv.DictReader(f):
        genes[row['gene-symbol']] = row
        ensembl[row['ensembl-id']] = row
assert len(genes) > 500, 'SFARI download is incomplete or has an unexpected format'
names = dict(line.split() for line in (REF / 'to_standard.tsv').read_text().splitlines())
reference_lengths = {f[0]: int(f[1]) for line in (REF/'GRCh38.fna.fai').read_text().splitlines()
                     if (f := line.split('\t'))}
fields = ['chromosome', 'position', 'ref', 'alt', 'gene', 'sfari_score', 'syndromic',
          'effect', 'predicted_impact', 'transcript', 'hgvs_c', 'hgvs_p',
          'genotype', 'depth', 'genotype_quality', 'allele_depths', 'alt_fraction',
          'filter', 'clinvar_significance', 'clinvar_condition', 'clinvar_review',
          'clinvar_allele_id', '1000G_AF', 'review_notes']
counts = {'all_variants': 0, 'pass_variants': 0, 'sfari_gene_variants': 0,
          'sfari_supported_autosomal_variants': 0, 'priority_review_variants': 0}
with (OUT / '571566.sfari_all.tsv').open('w') as all_file, \
     (OUT / '571566.sfari_priority_review.tsv').open('w') as priority_file:
    all_writer = csv.DictWriter(all_file, fields, delimiter='\t')
    priority_writer = csv.DictWriter(priority_file, fields, delimiter='\t')
    all_writer.writeheader()
    priority_writer.writeheader()
    with gzip.open(OUT / '571566.annotated.vcf.gz', 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            counts['all_variants'] += 1
            counts['pass_variants'] += c[6] == 'PASS'
            info = dict(x.split('=', 1) if '=' in x else (x, True)
                        for x in c[7].split(';'))
            fmt = dict(zip(c[8].split(':'), c[9].split(':')))
            gt = fmt.get('GT', '.')
            if not any(x not in {'0', '.'} for x in gt.replace('|', '/').split('/')):
                continue
            annotations = {}
            for annotation in str(info.get('ANN', '')).split(','):
                ann = annotation.split('|')
                if len(ann) < 11 or ann[0] != c[4]:
                    continue
                gene = genes.get(ann[3]) or ensembl.get(ann[4].split('.')[0])
                if gene:
                    symbol = gene['gene-symbol']
                    rank = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}
                    if symbol not in annotations or rank.get(ann[2], 9) < rank.get(annotations[symbol][1][2], 9):
                        annotations[symbol] = (gene, ann)
            if not annotations:
                continue
            counts['sfari_gene_variants'] += 1
            chrom = names.get(c[0], c[0])
            notes = []
            ad = fmt.get('AD', '.')
            fraction = None
            try:
                depths = list(map(int, ad.split(',')))
                fraction = depths[1] / sum(depths) if sum(depths) else None
            except (ValueError, IndexError):
                pass
            if chrom in {'X', 'Y'}:
                notes.append('sex-chromosome ploidy unconfirmed; genotype provisional')
            if c[6] != 'PASS':
                notes.append('fails variant or genotype support filters')
            alleles = gt.replace('|', '/').split('/')
            if '0' in alleles and '1' in alleles and (fraction is None or not 0.2 <= fraction <= 0.8):
                notes.append('heterozygous allele balance requires review')
            pop_af = str(info.get('POP_AF', '.'))
            try:
                common = max(float(x) for x in pop_af.split(',')) >= 0.01
            except ValueError:
                common = False
                notes.append('population frequency unavailable; does not establish rarity')
            if common:
                notes.append('frequency >=1% in this limited 1000 Genomes resource')
            if 'WARNING' in str(info.get('ANN', '')) or 'ERROR' in str(info.get('ANN', '')):
                notes.append('annotation warning/error; inspect ANN field')
            supported = not notes or all(x.startswith(('population frequency', 'frequency >=')) for x in notes)
            counts['sfari_supported_autosomal_variants'] += supported
            priority_variant = False
            for symbol, (gene, ann) in annotations.items():
                clnsig = str(info.get('CLNSIG', '.'))
                row = dict(zip(fields, [chrom, c[1], c[3], c[4], symbol,
                    gene['gene-score'], gene['syndromic'], ann[1], ann[2], ann[6],
                    ann[9], ann[10], gt, fmt.get('DP', '.'), fmt.get('GQ', '.'),
                    ad, '.' if fraction is None else f'{fraction:.3f}', c[6],
                    clnsig, info.get('CLNDN', '.'), info.get('CLNREVSTAT', '.'),
                    info.get('ALLELEID', '.'), pop_af, '; '.join(notes)]))
                all_writer.writerow(row)
                clinical_flag = 'pathogenic' in clnsig.lower()
                if supported and (clinical_flag or (ann[2] == 'HIGH' and not common)):
                    priority_writer.writerow(row)
                    priority_variant = True
            counts['priority_review_variants'] += priority_variant

coverage_lines = []
for depth in (1, 10, 20):
    path = OUT / f'qc/571566.coverage.min{depth}.tsv'
    with path.open() as f:
        data = list(csv.DictReader(f, delimiter='\t'))
    primary = [r for r in data if r['#rname'] in names]
    size = sum(reference_lengths[name] for name in names)
    covered = sum(int(r['covbases']) for r in primary)
    mean = sum(float(r['meandepth']) * (int(r['endpos']) - int(r['startpos']) + 1) for r in primary) / size
    coverage_lines.append(f'- Nuclear primary bases at >= {depth}x: {100*covered/size:.2f}%.')
    if depth == 1:
        coverage_lines.append(f'- Mean nuclear primary depth after quality/duplicate filtering: {mean:.2f}x.')
with (OUT/'qc/571566.sfari_gene_coverage.tsv').open() as f:
    coding = list(csv.DictReader(f, delimiter='\t'))
coverage_lines.append(f'- SFARI genes with MANE coding coverage measurements: {len(coding)}.')
coverage_lines.append(f'- Genes with >=90% of measured coding/flank bases at >=20x: '
                      f'{sum(float(r["percent_ge20"]) >= 90 for r in coding)}.')
(OUT / 'summary.json').write_text(json.dumps(counts, indent=2) + '\n')
(OUT / 'REPORT.md').write_text('''# Single-end WES research analysis: 571566

This is a candidate-variant screen, not an autism diagnosis or a clinical genetic report.
Gene membership and computational effect predictions do not establish pathogenicity.
The sample is WES: low whole-genome coverage is expected and is not an exome-quality
failure. Coding-region coverage is the relevant measure. The exact capture-kit BED
was unavailable; annotation-derived exon targets cannot establish kit-specific
capture completeness or on-target enrichment.

## Coverage

''' + '\n'.join(coverage_lines) + '''

Coverage excludes duplicate, secondary, supplementary, QC-failed and unmapped reads;
minimum mapping and base quality are both 20. Shallow or uncovered regions cannot
support negative findings. Single-end duplicate marking can also label independent
reads with identical alignment starts; optical duplicate counts are unavailable.

## Variant counts

''' + '\n'.join(f'- {k}: {v:,}' for k, v in counts.items()) + '''

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
''')
