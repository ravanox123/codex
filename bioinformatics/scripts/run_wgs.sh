#!/usr/bin/env bash
# Local, sequential, restartable single-end WGS analysis. Original BAM is retained.
set -euo pipefail
cd /home/ravanox/codex/bioinformatics
out="$PWD/results/wgs"
ref="$PWD/data/reference/GRCh38.fna"
sam="$PWD/tools/samtools-env/bin/samtools"
bcf="$PWD/tools/variants-env/bin/bcftools"
gatk="$PWD/tools/gatk-4.7.0.0/gatk"
snpeff="$PWD/tools/variants-env/bin/snpEff"
trap 'printf "FAILED at line %s (%s)\n" "$LINENO" "$(date -Is)" > "$out/STATUS.txt"' ERR
stage() { printf '%s %s\n' "$(date -Is)" "$*" | tee "$out/STATUS.txt"; }
run_gatk() { "$gatk" --java-options "-Xmx2g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=$out/tmp" "$@"; }
stage 'Waiting for duplicate marking and compatible reference resources'
while [[ ! -f "$out/markduplicates.complete" || ! -f "$out/resources.complete" ]]; do sleep 10; done
for source in 1000G_phase1.snps.high_confidence.hg38 Mills_and_1000G_gold_standard.indels.hg38 clinvar.GRCh38; do
    [[ -f "data/reference/$source.matched.vcf.gz.tbi" ]] || "$bcf" index -t "data/reference/$source.matched.vcf.gz"
done
if [[ ! -f "$out/qc/coverage.complete" ]]; then
    stage 'Measuring duplicate-filtered coverage (mapping/base quality >=20)'
    for depth in 1 10 20; do
        "$sam" coverage -q 20 -Q 20 --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY \
            --min-depth "$depth" -o "$out/qc/571566.coverage.min${depth}.tsv" "$out/571566.markdup.bam"
    done
    touch "$out/qc/coverage.complete"
fi
if [[ ! -f "$out/bqsr.complete" ]]; then
    stage 'Learning base-quality recalibration from known SNPs and indels'
    run_gatk BaseRecalibrator -R "$ref" -I "$out/571566.markdup.bam" \
        --default-platform UNKNOWN \
        --known-sites data/reference/1000G_phase1.snps.high_confidence.hg38.matched.vcf.gz \
        --known-sites data/reference/Mills_and_1000G_gold_standard.indels.hg38.matched.vcf.gz \
        -O "$out/571566.recalibration.table" > "$out/logs/base_recalibrator.log" 2>&1
    stage 'Applying base-quality recalibration'
    run_gatk ApplyBQSR -R "$ref" -I "$out/571566.markdup.bam" \
        --bqsr-recal-file "$out/571566.recalibration.table" -O "$out/571566.recal.bam" \
        > "$out/logs/apply_bqsr.log" 2>&1
    "$sam" quickcheck -v "$out/571566.recal.bam"
    touch "$out/bqsr.complete"
fi
# Nuclear primary chromosomes only: avoid treating alternate loci as independent calls.
# X/Y diploid genotypes are provisional until the sample's ploidy is supplied.
if [[ ! -f "$out/calling.complete" ]]; then
    stage 'Preparing WES exon targets (capture-kit BED unavailable)'
    while [[ ! -f "$out/qc/gene_coverage_full.complete" ]]; do sleep 10; done
    python3 scripts/build_wes_targets.py > "$out/logs/wes_targets.log"
    stage 'Calling WES exon-region SNPs and small indels (HaplotypeCaller)'
    run_gatk HaplotypeCaller -R "$ref" -I "$out/571566.recal.bam" \
        -L data/reference/wes.exons_padded100.bed --sample-ploidy 2 \
        --native-pair-hmm-threads 2 -O "$out/571566.raw.vcf.gz" \
        > "$out/logs/haplotypecaller.log" 2>&1
    touch "$out/calling.complete"
fi
if [[ ! -f "$out/filtering.complete" ]]; then
    stage 'Normalizing and filtering variant calls; retaining failed calls with labels'
    "$bcf" norm -f "$ref" -c e -m -any "$out/571566.raw.vcf.gz" -Oz -o "$out/571566.normalized.vcf.gz"
    "$bcf" filter -m + -s SNPQuality -e 'TYPE="snp" && (QUAL<30 || INFO/QD<2 || INFO/FS>60 || INFO/SOR>3 || INFO/MQ<40 || INFO/MQRankSum< -12.5 || INFO/ReadPosRankSum< -8)' \
        "$out/571566.normalized.vcf.gz" -Ou \
        | "$bcf" filter -m + -s IndelQuality -e 'TYPE="indel" && (QUAL<30 || INFO/QD<2 || INFO/FS>200 || INFO/ReadPosRankSum< -20)' -Ou \
        | "$bcf" filter -m + -s LowGenotypeSupport -e 'FMT/DP="." || FMT/GQ="." || FMT/DP<10 || FMT/GQ<20' \
            -Oz -o "$out/571566.filtered.vcf.gz"
    "$bcf" index -t "$out/571566.filtered.vcf.gz"
    touch "$out/filtering.complete"
fi
if [[ ! -f "$out/annotation.complete" ]]; then
    bash scripts/annotate_chunks.sh
fi
stage 'COMPLETE: small-variant research analysis; see REPORT.md for limitations and candidates'
