"""Build a local, unsigned consultation package; never send data or order tests."""
import csv,hashlib,html,json,re,shutil,zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/wgs/clinical_handoff'

TEXT='''# Clinical genetics consultation package — sample 571566

Research handoff draft, prepared 20 October 2026. Not a diagnosis, laboratory requisition, signed referral or completed clinical consultation.

## Status of the four requested steps

| Step | Prepared locally | Still needed |
|---|---|---|
| 1. Clinical genetics consultation | Genom Clinic selected; consultation request, report and variant table prepared | Patient details and appointment booking; Baku branch assumed from available context |
| 2. Confirm principal findings | Exact alleles and confirmation questions below | Clinician/laboratory acceptance, specimen collection and validated testing |
| 3. Clinical and family history | Structured form below | Patient/parent answers, records and clinician examination |
| 4. VWF review | Finding and assessment questions below | Bleeding history and clinician-directed assessment |

## Consultation request — draft to take to the clinic

To: Genom Clinic. The Baku clinic is assumed: [official website](https://genomclinic.az/), Mərkəzi Bulvar küç. 5, Ağ Şəhər, Bakı; +994 50 999 03 47 / +994 50 999 03 27; info@genomclinic.az. These are public contact details; the clinic has not been contacted and acceptance of external research data has not been verified.

Please assess sample 571566 following a research single-end whole-exome analysis for an autism-related investigation. The requester reports that clinical/family history is unavailable; this does not mean symptoms or family findings are absent. Age, pedigree and previous clinical testing have not been supplied. Two apparently homozygous findings in FDXR and WNT10A, and one heterozygous VWF finding, passed a conservative computational quality screen. Please review sample identity, independently confirm clinically relevant findings and evaluate phenotype/inheritance. These findings have not established a molecular explanation for autism.

This request is unsigned. It is prepared for the patient or caregiver to discuss with a clinician; no clinic has been contacted and no test has been ordered.

## Principal variants for laboratory discussion

| Gene | Transcript / protein | GRCh38 allele | Research call | Read support | Database evidence |
|---|---|---|---|---|---|
| FDXR | NM_024417.5:c.1156C>T; p.Arg386Trp | chr17:74863914 G>A; NC_000017.11:g.74863914G>A | Apparently homozygous; GQ 67 | Caller REF/ALT 1/35; local BQ20/MQ20 recount 1/34 | [ClinVar 520994](https://www.ncbi.nlm.nih.gov/clinvar/variation/520994/): pathogenic/likely pathogenic, multiple submitters |
| WNT10A | NM_025216.3:c.382C>T; p.Arg128Ter | chr2:218889989 C>T; NC_000002.12:g.218889989C>T | Apparently homozygous; GQ 99 | Caller REF/ALT 1/65; local recount 1/65 | [ClinVar 464181](https://www.ncbi.nlm.nih.gov/clinvar/variation/464181/): pathogenic/likely pathogenic in retained review |
| VWF | NM_000552.5:c.2561G>A; p.Arg854Gln | chr12:6034812 C>T; NC_000012.12:g.6034812C>T | Heterozygous; GQ 99 | Caller and local recount REF/ALT 48/39 | [ClinVar 296](https://www.ncbi.nlm.nih.gov/clinvar/variation/296/): expert-panel pathogenic for recessive type 2N disease in retained review |

The independent read counts are computational checks of the same sequencing specimen, not independent laboratory confirmation. ClinVar classifications concern variants, not a diagnosis in this patient. No independent ACMG/AMP classification has been completed. FDXR p.Arg392Trp is an alternative-transcript name for the same genomic allele.

## Confirmation questions for the receiving laboratory

| Question | Requested clarification |
|---|---|
| Sample identity | Confirm patient/sample linkage and laboratory specimen requirements; collection has not been arranged |
| Exact sequence calls | Confirm the listed alleles and genotypes through a clinically validated method selected by the laboratory |
| Apparently homozygous calls | Determine whether genotype interpretation is secure; consider allele dropout or a deletion/other copy-number explanation when appropriate |
| Inheritance | After confirmation, consider parental testing/segregation with clinical counseling and appropriate consent |
| Clinical classification | Review current variant evidence, gene–disease validity, phenotype and inheritance; do not substitute this draft for a laboratory interpretation |
| Remaining diagnostic scope | Consider gaps in exome coverage and unassessed variant classes as appropriate; this shortlist is not a complete clinical exome interpretation |

## Clinical information to complete

The requester reports no available history. All unprovided history items below remain unknown, not negative. Genom Clinic can collect history during consultation. For each item, record present, absent or unknown, with age of onset, dates and supporting records where available. This is a history form, not a symptom diagnosis checklist.

| History item | Answer / details |
|---|---|
| Patient identifier and relationship to sample 571566 | Not provided |
| Age / date of birth; sex recorded for clinical interpretation | Not provided |
| Reason for testing; established diagnoses and evaluating clinician | Not provided |
| Development: first words, walking, learning and speech concerns | Not provided |
| Developmental regression, seizures, movement or balance concerns | Not provided |
| Hearing concerns and previous audiology results | Not provided |
| Vision concerns and previous ophthalmology results | Not provided |
| Missing or unusually shaped teeth; previous dental imaging | Not provided |
| Skin, hair, nail or sweating concerns | Not provided |
| Easy bruising, nose/gum bleeding, prolonged bleeding after cuts | Not provided |
| Dental/surgical bleeding, heavy menstrual bleeding if applicable | Not provided |
| Prior anemia, transfusions or abnormal clotting results | Not provided |
| Current medicines and supplements relevant to bleeding assessment | Not provided |
| Similar developmental, hearing, vision, dental or bleeding findings in relatives | Not provided |
| Known parental relatedness, ancestry and affected relatives | Not provided |
| Previous genetic tests, laboratory names, dates and reports | Not provided |
| Availability of biological parents for clinical segregation discussion | Not provided |
| City/country, preferred clinic and referring clinician | Genom Clinic; Baku, Azerbaijan assumed; referring clinician not provided |

## VWF assessment for clinician review

The heterozygous p.Arg854Gln finding alone does not establish recessive type 2N von Willebrand disease. No second pathogenic coding allele was established; the research analysis does not exclude a second allele. Interpretation should incorporate bleeding and family history. [ClinVar evidence](https://www.ncbi.nlm.nih.gov/clinvar/variation/296/).

| Topic | Question for the clinician |
|---|---|
| Bleeding history | Is there a personal or family bleeding history warranting a hematology assessment? |
| Laboratory assessment if clinically indicated | Would VWF antigen, platelet-dependent VWF activity and factor VIII activity help assess the phenotype? |
| Suspected type 2N | Would specialist VWF–factor VIII binding testing or additional genetic evaluation be appropriate? These are clinician-selected options, not tests already ordered. |
| Interpretation | A heterozygous finding does not itself establish type 2N disease; heterozygous individuals may be asymptomatic, and symptoms require separate assessment. |

Assessment context: [ASH/ISTH/NHF/WFH diagnosis guideline](https://pmc.ncbi.nlm.nih.gov/articles/PMC7805340/). No treatment, medication or screening changes are prescribed by this package.

## Important limitations and attached evidence

| Area | Status |
|---|---|
| Autism/developmental explanation | Not established |
| Confirmation / phase / parental segregation | Not completed |
| Copy-number, structural, repeat-expansion and mitochondrial assessment | Not performed in this review |
| SYNGAP1 coding/flank bases reaching 20× | 0.17%; substantial gap |
| SHANK3 coding/flank bases reaching 20× | 87.62% in earlier supplementary measurement |
| Remaining variant interpretation | 128 other MANE HIGH candidates remain an unclassified research queue |
| CNN2 substitutions | Fail stricter quality checks; excluded from principal findings |
| Files transmitted externally | None by this package builder |

The attached report, quality audit summary and eight-site variant table preserve the research evidence. Visual-review statuses remain as recorded in the export; an IGV screenshot is not equivalent to a signed manual review. Raw BAM/FASTQ and full VCF files are not included in this consultation archive. They remain available locally if a receiving laboratory requests them through an agreed transfer method.
'''

def render(text):
    def inline(s):
        s=html.escape(s)
        return re.sub(r'\[([^]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',s)
    parts=[];table=False;rows=0
    for line in text.splitlines():
        if line.startswith('|'):
            cells=[c.strip() for c in line.strip('|').split('|')]
            if all(re.fullmatch('[: -]+',c) for c in cells):continue
            if not table:parts.append('<table>');table=True;rows=0
            tag='th' if rows==0 else 'td';parts.append('<tr>'+''.join(f'<{tag}>{inline(c)}</{tag}>' for c in cells)+'</tr>');rows+=1;continue
        if table:parts.append('</table>');table=False
        if not line.strip():continue
        m=re.match(r'^(#{1,6}) (.*)',line)
        parts.append(f'<h{len(m[1])}>{inline(m[2])}</h{len(m[1])}>' if m else '<p>'+inline(line)+'</p>')
    if table:parts.append('</table>')
    return '<!doctype html><html lang="en"><meta charset="utf-8"><title>571566 — Clinical handoff</title><style>body{font:14px/1.5 Arial;margin:30px;color:#18373e}h1,h2{color:#125661}table{border-collapse:collapse;width:100%;font-size:12px}th,td{border:1px solid #c8d8d6;padding:9px;vertical-align:top;text-align:left}th{background:#e8f2f0}tr:nth-child(even){background:#f7f9f8}h2{margin-top:28px}button{padding:10px}@page{size:A4 landscape;margin:12mm}@media print{button{display:none}tr{break-inside:avoid}h2{break-after:avoid}}</style><button onclick="window.print()">Print / Save as PDF</button>'+''.join(parts)+'</html>'

def main():
    OUT.mkdir(exist_ok=True)
    (OUT/'START_HERE.md').write_text(TEXT)
    (OUT/'START_HERE.html').write_text(render(TEXT))
    sources={
        'CLINICAL_STYLE_REPORT.md':ROOT/'results/wgs/exome_review/CLINICAL_STYLE_REPORT.md',
        'CLINICAL_STYLE_REPORT.html':ROOT/'results/wgs/exome_review/CLINICAL_STYLE_REPORT.html',
        'REVIEW.md':ROOT/'results/wgs/exome_review/REVIEW.md',
        'CLINICAL_HANDOFF.md':ROOT/'results/wgs/exome_review/CLINICAL_HANDOFF.md',
        '571566_variant_review.csv':ROOT/'results/wgs/igv_review/571566_variant_review.csv',
        '571566_variant_review.tsv':ROOT/'results/wgs/igv_review/571566_variant_review.tsv',
        '571566_variant_review.html':ROOT/'results/wgs/igv_review/571566_variant_review.html',
        'QUALITY_REPORT.md':ROOT/'results/wgs/quality_review/QUALITY_REPORT.md',
        'quality_summary.json':ROOT/'results/wgs/quality_review/summary.json',
        'report_findings_quality.tsv':ROOT/'results/wgs/quality_review/report_findings_quality.tsv',
        'assertion_SNV_read_support.json':ROOT/'results/wgs/exome_review/assertion_SNV_read_support.json',
    }
    for name,source in sources.items():shutil.copy2(source,OUT/name)
    history=TEXT.split('## Clinical information to complete',1)[1].split('## VWF assessment',1)[0]
    with (OUT/'clinical_history_form.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['History item','Answer / details'])
        for line in history.splitlines():
            if line.startswith('|') and not line.startswith(('|---','| History item')):writer.writerow([c.strip() for c in line.strip('|').split('|')])
    files=[OUT/n for n in ['START_HERE.md','START_HERE.html','clinical_history_form.csv',*sources]]
    manifest={'sample':'571566','status':'unsigned research handoff; no appointment or test ordered','preferred_clinic':'Genom Clinic; Baku branch assumed','clinical_history':'requester reports history unavailable; no negative findings inferred','files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    archive=OUT.parent/'571566_clinical_handoff.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in [*files,OUT/'manifest.json']:z.write(p,p.name)
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        assert len(z.namelist())==len(files)+1
    assert (OUT/'START_HERE.html').read_text().count('<table>')==6
    print(f'Created {archive} ({archive.stat().st_size:,} bytes); {len(files)+1} files verified.')

if __name__=='__main__':main()
