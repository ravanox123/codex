# WES Lab

A native Python desktop application for local WES research: searchable variants, SFARI gene overlap, transcript context, interactive graphs, measured coverage, saved projects and exports. It does not calculate an autism probability or make a diagnosis.

## Start

A dedicated Python 3.12 environment with Tk and proper desktop font rendering is installed at `tools/wes-lab-env`:

```bash
./scripts/start_wes_lab.sh
```

Open **desktop-shortcuts/WES-Lab.desktop**, or load the existing sample directly:

```bash
./scripts/start_wes_lab.sh --example
```

To reopen a particular project:

```bash
./scripts/start_wes_lab.sh --project results/wes_lab/571566
```

On another computer, Python 3.10+ with Tk 8.6 is required. There are no third-party Python dependencies for the viewer. Run `python -m wes_lab`, or install this package and run `wes-lab`. The launcher uses the dedicated app environment because this computer's system Python does not include Tk. The annotation environments are separate.

## Add WES data

1. Click **Add WES / VCF** and select a `.vcf` or `.vcf.gz` file.
2. Choose the sample if the VCF contains several samples. Choose the actual genome build; keep Unknown if it has not been established.
3. Import runs on a background thread into a SQLite project under `results/wes_lab/`. The source file is not copied or edited. A cancelled or failed import does not replace a previously completed database.
4. Use **Variants** to search genes, HGVS or chromosome:position; filter by PASS, SFARI, representative effect, chromosome and review group. Double-click a row for all transcript annotations.
5. Attach a measured gene coverage TSV in **Coverage**. Required columns are `gene`, `mean_depth` and `percent_ge20`; optional `coding_plus_flank_bases` supplies interval lengths. A new attachment replaces the previous coverage table. This file must describe the same sample; the app cannot verify specimen identity from a coverage table alone.

SnpEff `ANN` and VEP `CSQ` are supported. An unannotated VCF can be imported for call-quality graphs, but gene-based review remains unavailable until annotations are provided. Only selected-sample nonreference alleles are imported when a complete genotype is available; unknown/partial genotypes are retained with an explicit support warning. Multi-ALT genotypes, AD and population frequencies are matched to the individual ALT allele. Original genotypes are preserved.

Bundled SFARI and GRCh38 MANE transcript resources are offline snapshots. The SFARI CSV was downloaded on 2026-10-16; provenance and source hash are in `wes_lab/assets/resources.json`. The app does not update public databases automatically. GRCh38 MANE ID comparisons are used only when GRCh38 is selected, while explicit VEP MANE Select tags can be used on other builds. Transcript versions may differ; a matching ID is not clinical transcript validation.

Representative effects prefer a matched MANE transcript before the strongest predicted effect. Alternative transcript effects remain visible. “Conflicting classifications of pathogenicity” is a separate group, never automatically treated as a pathogenic assertion. Population `INFO/POP_AF` is used when supplied; generic `INFO/AF` is not assumed to be a population frequency because callers often use it for the sample callset. Missing frequency never establishes rarity. Record-level ClinVar assertions on unsplit multi-ALT records require allele matching and are flagged for support review.

## Graphs and exports

**Overview** includes clickable transcript-effect and chromosome charts, variant types, variant-site depth, clinical assertions and measured gene-coverage bins. Hover over bars for counts; use **Save graph** for standalone SVG files. Charts summarize the selected sample's imported alleles, including filtered calls, unless their captions say otherwise. Variant-site DP is not an exome coverage measurement.

**Reports** exports a standalone HTML report with inline SVG charts; use a browser's Print function for PDF. CSV exports stream all matching variants rather than only the current page. Pagination keeps the table responsive. Coverage can be exported separately as TSV.

Reports and exports contain genomic information and remain on local disk. The app makes no network requests while importing or reviewing files. **Search public ClinVar** deliberately opens a gene/HGVS search in a browser when clicked. No BAM, FASTQ or VCF upload is performed.

## Local BAM and annotation workflow

**Run analysis** supports:

- Calling SNPs/small indels from an **analysis-ready BAM** (coordinate sorted, indexed, read groups for one sample, compatible with the installed reference).
- Annotating a single-sample, nuclear-primary, accession-named GRCh38 VCF with the installed resource bundle.

The workflow uses local GATK, samtools, bcftools, SnpEff GRCh38 Ensembl 115, ClinVar and the limited frequency resource under `tools/` and `data/`. It requires all resources reported by `wes_lab.pipeline.probe_bundle()`; they are already present here. This is not a portable downloader or a general pipeline for arbitrary reference assemblies. The viewer itself can import VCFs without those tools.

BAM calling is sequential by chromosome, intersected with the installed Ensembl exon+100bp targets, with 100bp assembly padding and trimming back to target start positions. The capture-kit BED is not available, so the supplied targets are approximate. Calling uses a 2 GB Java heap and one native caller thread. SNP, indel and genotype hard filters follow the existing research workflow. X/Y are called with provisional diploid ploidy and remain flagged.

Annotation runs in occupied 10 Mb sections, one at a time, retaining chromosome-specific transcripts with a 3 GB Java heap, serial GC and one JVM processor. Full transcript context on that chromosome is available. Completion checkpoints and configuration/resource fingerprints prevent accidental reuse of an incompatible workflow. Final merging checks the complete chromosome/position/REF/ALT multiset, including duplicate multiplicities. The files are published only after that check passes. A workflow lock prevents another runner using the same folder. Stop terminates the running subprocess group; restart with the same input and folder to resume.

Do not start another memory-heavy annotation process while this workflow runs. In particular, WES Lab cannot detect a separate manually launched pipeline using another folder.

**Not included:** FASTQ alignment, read preparation/duplicate marking/BQSR, CNV/SV calling, repeat expansions, mitochondrial interpretation, mosaicism analysis, parental phasing or a clinically validated classifier. Input BAM preparation and recalibration suitability must already have been assessed. The workflow's biological accuracy has not been validated by the software tests; no new full GATK run was performed as part of building the viewer.

## Existing sample 571566

The saved viewer project is `results/wes_lab/571566`. It contains the complete annotated input and measured MANE coverage, supplemented with the existing SFARI/Ensembl coverage where MANE genes were absent. **Reports** also opens the original 18-candidate review and the broader exome review/handoff. Those manual reports distinguish FDXR/WNT10A calls requiring clinical confirmation from weaker or artifact-prone findings; the viewer's automatic flags do not replace those reviews.

## Checks

```bash
python3 -m unittest discover -s tests -v
```

Tests cover allele-specific multi-sample/multi-ALT parsing, VEP matching, MANE preference, clinical-conflict handling, missing frequencies, cancelled/failed imports, coverage validation, HTML escaping and subprocess cancellation. The real 97,174-record sample is also imported and cross-checked against the original callset. It contains 110 spanning-deletion `*` placeholders, which are excluded from standalone variant displays, leaving 97,064 displayed alleles and 73,290 PASS alleles. Placeholder exclusions are recorded in project metadata.

## Person projects and file downloads

**Projects** now treats each person as one project. Create a project with a name
or anonymous identifier and optional notes. Each person can hold multiple WES
analyses; the **Files → Analyses** tab switches between them. Names are labels,
not verification that different files belong to the same biological person.

Use **Files** to:

- **Download links:** paste direct HTTP/HTTPS links, one per line. Paired FASTQ
  URLs can be queued together. For a single link, optionally supply a filename
  and expected SHA-256 checksum.
- **Link existing files:** register files in place without duplicating large
  FASTQ/BAM files. Their original location must remain accessible.
- **Copy files here:** stream copies into this person's folder; original files
  remain in place. Name collisions are rejected rather than overwritten.
- **Use selected file in analysis:** import a VCF into a new analysis under this
  person, or select a BAM as the input to the existing prepared-BAM workflow.
- **Open folder:** inspect the organized files in the file manager.

Person folders live in `results/wes_lab/people/<identifier>_<unique-id>/`:

```text
person.json          person label, notes, file registry and analysis links
files/fastq/         managed FASTQ/FQ files, including .gz
files/bam/           BAM/CRAM files and BAI/CRAI indexes
files/vcf/           VCF/BCF files and TBI/CSI indexes
files/other/         other supplied files
analyses/            imported variant databases and workflow outputs
reports/             default destination for report/table exports
downloads/           persistent download jobs and unpublished .part files
```

The **Downloads** tab shows queued, paused, failed and completed transfers with
sizes and source locations. Select a row for error details. **Stop task** pauses
the current transfer; **Resume selected download** or **Resume all pending**
continues saved jobs after an interruption or app restart. Resume appends bytes
only when the server supplies an ETag or Last-Modified validator and a valid
HTTP range response; otherwise the partial transfer restarts. Interrupted
copies are removed; local copy operations are not resumable.

Downloads stream to partial files and appear in Files only after completion
and SHA-256 calculation. If an expected checksum is provided, it must match
before publication. A computed checksum without an independently supplied
expected value identifies the file but does not authenticate its contents.
FASTQ/BAM format validity and biological sample identity are not automatically
verified by downloading a file. Downloading FASTQ is supported; FASTQ alignment
is still outside the analysis workflow.

Existing file contents are never silently overwritten. Downloads accept direct
HTTP/HTTPS file URLs rather than portal/login pages, FTP links or SRA accession
numbers. Signed HTTPS links work; their query strings are kept only in private
local resume-job metadata and omitted from visible source URLs. No sequencing
files are uploaded. Use local authorization and access-controlled storage for
sensitive datasets as appropriate.

The existing 571566 analysis is registered under a person project with links to
the existing BAM, VCF and reviews; its large source files are not moved or copied.
Old standalone analysis folders can still be opened and are registered under a
person automatically. New imports and workflows default to the selected person's
analysis directory; report exports default to that person's reports directory.

The test suite now includes person isolation, file copying/linking, path escape
protection, and real local HTTP transfers with pause/resume, checksum mismatch,
invalid ranges, truncated responses and filename collisions. Those network tests
bind a temporary localhost server and transfer only synthetic data.
# IGV read review for sample 571566

Run `scripts/start_variant_review.sh`, or open the existing 571566 analysis in WES Lab and choose **Variants → Inspect reads / IGV**. Select a row and choose **Open selected in IGV**, or double-click it. The control window searches all quality-audited records and initially shows the eight curated review sites. Other variants have no curated pathogenicity label in this controller.

IGV Desktop 2.19.8 displays the local GRCh38 FASTA, recalibrated BAM, retained/review VCF tracks and report-site intervals. Accession chromosome names are mapped to the exact local reference. IGV uses its bundled graphical Java runtime and a 2 GB heap instead of its default 8 GB; only the selected region is requested. No sequencing files are uploaded. The controller is specifically configured for 571566 and refuses to attach these reads to another WES Lab analysis.

Visual review statuses and notes are saved under `results/wgs/igv_review/visual_review.json`; **Export reviews TSV** writes a table in that folder. These are review notes, not new clinical classifications. The session XML and IGV log are in the same directory. Close IGV separately when finished; closing the control window leaves IGV open. If no matching variant is selected, choose a row first. **All variants**, **Quality retained** and **Needs quality review** show up to 500 matches; narrow the search to inspect further records.

The quality-screen provenance and thresholds are documented in `results/wgs/quality_review/QUALITY_REPORT.md`. Curated pathogenicity labels come from the retained review and are not automatically refreshed from ClinVar. Official IGV download: https://igv.org/doc/desktop/DownloadPage/.

