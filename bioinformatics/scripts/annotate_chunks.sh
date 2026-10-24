#!/usr/bin/env bash
# Sequential, restartable annotation. Windows assign variants by POS, preserving coordinates.
set -euo pipefail
cd /home/ravanox/codex/bioinformatics
out="$PWD/results/wgs"
chunks="$out/annotation_chunks"
bcf="$PWD/tools/variants-env/bin/bcftools"
mkdir -p "$chunks"
exec 8>"$out/annotation_chunks.lock"
flock -n 8 || { echo 'Chunk annotation is already running' >&2; exit 1; }
trap 'printf "FAILED: chunk annotation at line %s (%s); see logs/annotation_chunks.log\n" "$LINENO" "$(date -Is)" > "$out/STATUS.txt"' ERR
# Fingerprint inputs and implementation so stale checkpoints cannot be reused.
sha256sum scripts/annotate_chunks.sh "$out/571566.filtered.vcf.gz" \
 data/reference/clinvar.GRCh38.matched.vcf.gz data/reference/1000G.af_available.vcf.gz \
 data/reference/to_standard.tsv data/reference/to_accessions.tsv data/reference/ensembl115.bed \
 data/snpeff/GRCh38.115/snpEffectPredictor.bin > "$chunks/inputs.current"
if [[ -f "$chunks/inputs.sha256" ]]; then
    cmp "$chunks/inputs.sha256" "$chunks/inputs.current" || { echo 'Inputs changed; use a fresh annotation_chunks directory.' >&2; exit 1; }
else
    mv "$chunks/inputs.current" "$chunks/inputs.sha256"
fi
python3 - <<'PY' > "$chunks/regions.tsv"
import subprocess
bcf='tools/variants-env/bin/bcftools'
vcf='results/wgs/571566.filtered.vcf.gz'
counts={line.split('\t')[0]:int(line.split('\t')[2]) for line in subprocess.check_output([bcf,'index','-s',vcf],text=True).splitlines()}
for line in open('data/reference/GRCh38.fna.fai'):
    chrom,length,*_=line.split()
    if not counts.get(chrom): continue
    for start in range(1,int(length)+1,10_000_000):
        print(chrom,start,min(start+9_999_999,int(length)),sep='\t')
PY
# Keep all transcripts on each current chromosome before SnpEff builds its forest.
python3 - <<'PYTRANSCRIPTS'
from pathlib import Path
from collections import defaultdict
out=Path('results/wgs/annotation_chunks')
names=dict(line.split() for line in Path('data/reference/to_standard.tsv').read_text().splitlines())
transcripts=defaultdict(set)
with open('data/reference/ensembl115.bed') as f:
    for line in f:
        chrom, start, end, feature = line.rstrip().split('\t')[:4]
        if feature.startswith('Transcript_'):
            transcripts[chrom].add(feature[len('Transcript_'):])
for accession, chrom in names.items():
    if transcripts[chrom]:
        (out / (accession+'.transcripts.txt')).write_text(''.join(t+'\n' for t in sorted(transcripts[chrom])))
PYTRANSCRIPTS
: > "$chunks/concat.list"
while IFS=$'\t' read -r chrom start end; do
    name="$chunks/${chrom}_${start}_${end}"
    if [[ ! -f "$name.complete" ]]; then
        printf 'Annotating %s:%s-%s (%s)\n' "$chrom" "$start" "$end" "$(date -Is)" | tee "$out/STATUS.txt"
        "$bcf" view -r "$chrom:$start-$end" --regions-overlap 0 "$out/571566.filtered.vcf.gz" -Oz -o "$name.input.vcf.gz"
        "$bcf" index -f -t "$name.input.vcf.gz"
        count=$("$bcf" index -n "$name.input.vcf.gz")
        if [[ "$count" != 0 ]]; then
            "$bcf" annotate -a data/reference/clinvar.GRCh38.matched.vcf.gz \
                -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT,INFO/ALLELEID "$name.input.vcf.gz" -Oz -o "$name.clinvar.vcf.gz"
            "$bcf" index -f -t "$name.clinvar.vcf.gz"
            "$bcf" annotate -a data/reference/1000G.af_available.vcf.gz -c INFO/POP_AF:=INFO/AF "$name.clinvar.vcf.gz" -Ou \
            | "$bcf" annotate --rename-chrs data/reference/to_standard.tsv -Ov \
            | tools/variants-env/bin/snpEff -Xmx3g -XX:+UseSerialGC -XX:ActiveProcessorCount=1 \
                -onlyTr "$chunks/$chrom.transcripts.txt" -dataDir "$PWD/data/snpeff" -noStats -noLog -nodownload GRCh38.115 2> "$name.snpeff.log" \
            | "$bcf" annotate --rename-chrs data/reference/to_accessions.tsv -Oz -o "$name.partial.vcf.gz"
            "$bcf" index -f -t "$name.partial.vcf.gz"
            [[ $("$bcf" index -n "$name.partial.vcf.gz") == "$count" ]]
            mv "$name.partial.vcf.gz" "$name.vcf.gz"
            mv "$name.partial.vcf.gz.tbi" "$name.vcf.gz.tbi"
        fi
        printf '%s\n' "$count" > "$name.complete"
    fi
    if [[ $(cat "$name.complete") != 0 ]]; then printf '%s\n' "$name.vcf.gz" >> "$chunks/concat.list"; fi
done < "$chunks/regions.tsv"
"$bcf" concat -f "$chunks/concat.list" -Oz -o "$out/571566.annotated.partial.vcf.gz"
"$bcf" index -f -t "$out/571566.annotated.partial.vcf.gz"
[[ $("$bcf" index -n "$out/571566.annotated.partial.vcf.gz") == $("$bcf" index -n "$out/571566.filtered.vcf.gz") ]]
# Ensure identical variant multisets, including duplicate multiplicities.
# Annotation may reorder different alleles sharing the same position.
cmp <("$bcf" query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$out/571566.filtered.vcf.gz" | LC_ALL=C sort) \
     <("$bcf" query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$out/571566.annotated.partial.vcf.gz" | LC_ALL=C sort)
mv "$out/571566.annotated.partial.vcf.gz" "$out/571566.annotated.vcf.gz"
mv "$out/571566.annotated.partial.vcf.gz.tbi" "$out/571566.annotated.vcf.gz.tbi"
"$bcf" view -f PASS "$out/571566.annotated.vcf.gz" -Oz -o "$out/571566.pass.vcf.gz"
"$bcf" index -f -t "$out/571566.pass.vcf.gz"
"$bcf" stats "$out/571566.annotated.vcf.gz" > "$out/qc/571566.vcf_stats.txt"
python3 scripts/report_wgs.py
touch "$out/annotation.complete"
printf 'COMPLETE: chunk annotation and report (%s)\n' "$(date -Is)" > "$out/STATUS.txt"
