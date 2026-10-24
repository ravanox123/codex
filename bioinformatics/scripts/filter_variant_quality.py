"""Auditable, conservative single-sample small-variant quality screen.

Original records are preserved in separate accepted/review VCFs. This is not
a calibrated clinical sensitivity/specificity model or pathogenicity filter.
"""
import collections
import csv
import gzip
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/wgs/quality_review'


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def assess(v, read_check=None):
    info = dict(x.split('=', 1) if '=' in x else (x, True) for x in v[7].split(';'))
    fmt = dict(zip(v[8].split(':'), v[9].split(':')))
    reasons, warnings = [], []
    if v[6] != 'PASS':
        reasons.append('original_filter:' + v[6])
    alt = v[4]
    if ',' in alt or not re.fullmatch('[ACGT]+', v[3]) or not re.fullmatch('[ACGT]+', alt):
        reasons.append('unsupported_or_placeholder_allele')
    gt = fmt.get('GT', '.').replace('|', '/').split('/')
    if any(a not in ('0', '1') for a in gt) or '1' not in gt or len(gt) not in (1, 2):
        reasons.append('missing_or_unsupported_genotype')
    def minimum(label, value, threshold):
        n = number(value)
        if n is None:
            reasons.append('missing_' + label)
        elif n < threshold:
            reasons.append(label + '_below_' + str(threshold))
    minimum('QUAL', v[5], 30)
    minimum('DP', fmt.get('DP'), 20)
    minimum('GQ', fmt.get('GQ'), 30)
    minimum('QD', info.get('QD'), 2)
    minimum('MQ', info.get('MQ'), 40)
    snp = len(v[3]) == len(alt) == 1
    for label, limit in [('FS', 60 if snp else 200), *([('SOR', 3)] if snp else [])]:
        n = number(info.get(label))
        if n is None:
            reasons.append('missing_' + label)
        elif n > limit:
            reasons.append(label + '_above_' + str(limit))
    # Rank sums can legitimately be unavailable for homozygous calls.
    for label, limit in [('ReadPosRankSum', -8 if snp else -20), *([('MQRankSum', -12.5)] if snp else [])]:
        n = number(info.get(label))
        if n is None:
            warnings.append(label + '_unavailable')
        elif n < limit:
            reasons.append(label + '_below_' + str(limit))
    balance = None
    try:
        ad = [int(a) for a in fmt.get('AD', '.').split(',')]
        if len(ad) != 2 or min(ad) < 0 or sum(ad) == 0:
            raise ValueError()
        balance = ad[1] / sum(ad)
        if ad[1] < 5:
            reasons.append('ALT_reads_below_5')
        if sorted(gt) == ['0', '1'] and not .25 <= balance <= .75:
            reasons.append('heterozygous_balance_outside_0.25_0.75')
        elif gt and set(gt) == {'1'} and balance < .9:
            reasons.append('homozygous_or_haploid_ALT_fraction_below_0.90')
    except ValueError:
        reasons.append('missing_or_invalid_AD')
    if read_check and read_check.get('alt_BQ20_MQ40', 0) < 5:
        reasons.append('independent_ALT_MQ40_reads_below_5')
    if read_check is None:
        warnings.append('no_independent_read_recount')
    return reasons, warnings, info, fmt, balance


def main():
    OUT.mkdir(exist_ok=True)
    names = dict(l.split() for l in (ROOT / 'data/reference/to_standard.tsv').read_text().splitlines())
    checks = json.loads((ROOT / 'results/wgs/exome_review/assertion_SNV_read_support.json').read_text())
    checks = {(str(r['chromosome']).removeprefix('chr'), str(r['position']), r['ref'], r['alt']): r for r in checks}
    counts, failures, lead_rows = collections.Counter(), collections.Counter(), []
    fields = ['chromosome', 'position', 'ref', 'alt', 'genes', 'GT', 'DP', 'GQ', 'AD', 'ALT_fraction', 'QUAL', 'QD', 'MQ', 'FS', 'SOR', 'original_filter', 'quality_status', 'failure_reasons', 'warnings']
    with gzip.open(ROOT / 'results/wgs/571566.annotated.vcf.gz', 'rt') as src, (OUT / 'accepted.vcf').open('w') as keep, (OUT / 'needs_review.vcf').open('w') as review, (OUT / 'variant_quality.tsv').open('w') as audit:
        writer = csv.DictWriter(audit, fields, delimiter='\t')
        writer.writeheader()
        for line in src:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    if len(line.rstrip().split('\t')) != 10:
                        raise ValueError('This screen requires exactly one sample')
                    extra = '##quality_screen="Research screen; see QUALITY_REPORT.md and variant_quality.tsv; original FILTER fields preserved"\n'
                    keep.write(extra); review.write(extra)
                keep.write(line); review.write(line)
                continue
            v = line.rstrip('\n').split('\t')
            chrom = names.get(v[0], v[0]).removeprefix('chr')
            check = checks.get((chrom, v[1], v[3], v[4]))
            reasons, warnings, info, fmt, balance = assess(v, check)
            if chrom in ('X', 'Y'):
                warnings.append('sex_chromosome_ploidy_not_validated')
            genes = sorted({a.split('|')[3] for a in info.get('ANN', '').split(',') if len(a.split('|')) > 3 and a.split('|')[0] == v[4]})
            row = dict(chromosome=chrom, position=v[1], ref=v[3], alt=v[4], genes=';'.join(genes), GT=fmt.get('GT', '.'), DP=fmt.get('DP', '.'), GQ=fmt.get('GQ', '.'), AD=fmt.get('AD', '.'), ALT_fraction='.' if balance is None else f'{balance:.5f}', QUAL=v[5], QD=info.get('QD', '.'), MQ=info.get('MQ', '.'), FS=info.get('FS', '.'), SOR=info.get('SOR', '.'), original_filter=v[6], quality_status='NEEDS_REVIEW' if reasons else 'ACCEPTED', failure_reasons=';'.join(reasons), warnings=';'.join(warnings))
            writer.writerow(row)
            (review if reasons else keep).write(line)
            counts['total'] += 1
            counts[row['quality_status']] += 1
            counts['original_PASS'] += v[6] == 'PASS'
            counts['original_PASS_now_review'] += v[6] == 'PASS' and bool(reasons)
            failures.update(reasons)
            if set(genes) & {'FDXR', 'WNT10A', 'VWF', 'CNN2', 'DNAH14', 'TRMT5', 'DOLK'} and ('Pathogenic' in info.get('CLNSIG', '') or any(a.split('|')[2:3] == ['HIGH'] for a in info.get('ANN', '').split(','))):
                lead_rows.append(row)
    with (OUT / 'report_findings_quality.tsv').open('w') as f:
        writer = csv.DictWriter(f, fields, delimiter='\t'); writer.writeheader(); writer.writerows(lead_rows)
    summary = {'counts': dict(counts), 'failure_reason_counts_nonexclusive': dict(failures), 'independent_read_check_loci': len(checks)}
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    report = ['# Variant quality screen — sample 571566', '', 'Research quality selection, not a clinical validation or pathogenicity classification. Original VCFs are unchanged. Rejected records are retained for review, not declared false positives.', '', '| Result | Records |', '|---|---:|']
    report += [f'| {k} | {v:,} |' for k, v in counts.items()]
    report += ['', '## Rules', '', '| Metric | Required |', '|---|---|', '| Original FILTER | PASS |', '| QUAL / sample DP / sample GQ | ≥30 / ≥20 / ≥30 |', '| QD / MQ | ≥2 / ≥40 (MQ extended to indels as a conservative local rule) |', '| FS | ≤60 SNP; ≤200 other small alleles |', '| SOR | ≤3 SNP; no extra indel SOR cutoff |', '| ReadPosRankSum | ≥−8 SNP; ≥−20 indel, when available |', '| MQRankSum | ≥−12.5 SNP, when available |', '| Allelic depths | Two nonnegative AD entries, at least 5 ALT reads |', '| Heterozygous ALT fraction | 0.25–0.75 |', '| Homozygous/haploid ALT fraction | ≥0.90 |', '| Existing independent SNV recount | At least 5 ALT reads at BQ20/MQ40 if a recount is available |', '| Missing values | Required metrics fail; rank sums may be unavailable and are flagged |', '', 'These are conservative research thresholds. DP, GQ, allele balance and the indel MQ rule are local choices, not a clinically calibrated model. Low-depth genuine variants, mosaics and imbalanced calls can be excluded. No maximum-depth filter was imposed on exome data. Population frequency, predicted impact and ClinVar significance do not affect quality acceptance. Sex-chromosome ploidy is unresolved and flagged. Multi-ALT, symbolic and spanning-deletion placeholder records are retained for review. The existing independent recount covers only a small subset; accepted variants have not all had read-level alignment validation.', '', 'Site-filter basis: [GATK hard-filtering guidance](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants).', '', '## Failure reasons (overlap possible)', '', '| Reason | Records |', '|---|---:|']
    report += [f'| {k} | {v:,} |' for k,v in failures.most_common()]
    report += ['', '## Findings from report genes', '', '| Gene | GRCh38 allele | GT | DP | GQ | ALT fraction | Result | Reasons |', '|---|---|---|---:|---:|---:|---|---|']
    report += [f"| {r['genes']} | {r['chromosome']}:{r['position']} {r['ref']}>{r['alt']} | {r['GT']} | {r['DP']} | {r['GQ']} | {r['ALT_fraction']} | {r['quality_status']} | {r['failure_reasons'] or '—'} |" for r in lead_rows]
    report += ['', 'Full per-record metrics, decisions and warnings: `variant_quality.tsv`. VCF subset FILTER columns retain their original meaning; the new screen decision is recorded in the TSV and by subset membership. The earlier clinical-style report remains a historical research draft; consult this quality addendum alongside it.']
    (OUT / 'QUALITY_REPORT.md').write_text('\n'.join(report) + '\n')
    assert counts['total'] == counts['ACCEPTED'] + counts['NEEDS_REVIEW']
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
