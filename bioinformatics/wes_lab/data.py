"""Streaming VCF import and bounded SQLite queries. No GUI or external dependencies."""
from __future__ import annotations
import csv,gzip,hashlib,html,json,math,os,re,sqlite3,threading,time,uuid
from contextlib import contextmanager
from collections import Counter
from pathlib import Path
from typing import Callable

ASSETS=Path(__file__).with_name('assets')
IMPACT_ORDER={'HIGH':0,'MODERATE':1,'LOW':2,'MODIFIER':3,'UNKNOWN':4}
class Cancelled(Exception):pass

def open_text(path):
    path=Path(path)
    return gzip.open(path,'rt',encoding='utf-8') if path.suffix.lower()=='.gz' else path.open(encoding='utf-8-sig')

def numeric(value):
    try:
        number=float(value)
        return number if math.isfinite(number) else None
    except (TypeError,ValueError):return None

def split_info(text):
    return dict(x.split('=',1) if '=' in x else (x,True) for x in text.split(';') if x and x!='.')

def classify_assertion(value):
    """Conflicting pathogenicity labels must never be mistaken for pathogenic assertions."""
    text=str(value or '')
    if 'conflicting' in text.lower():return 'Conflicting'
    tokens={x.strip().lower().replace(' ','_') for x in re.split(r'[/|,;]',text)}
    if tokens&{'pathogenic','likely_pathogenic'}:return 'Pathogenic / likely pathogenic'
    if tokens&{'benign','likely_benign'}:return 'Benign / likely benign'
    if tokens&{'uncertain_significance','variant_of_uncertain_significance'}:return 'Uncertain significance'
    if text in ('','.','None'):return 'Not annotated'
    return 'Other assertion'

def read_header(path):
    samples=[];csq=[];info_numbers={};reference='Not declared'
    with open_text(path) as f:
        for line in f:
            if line.startswith('##reference='):reference=line.strip().partition('=')[2]
            if line.startswith('##INFO=<'):
                match=re.search(r'ID=([^,>]+),Number=([^,>]+)',line)
                if match:info_numbers[match[1]]=match[2]
                if 'ID=CSQ,' in line and 'Format:' in line:
                    csq=line.split('Format:',1)[1].strip().rstrip('>"\r\n ').split('|')
            if line.startswith('#CHROM\t'):
                samples=line.rstrip('\r\n').split('\t')[9:]
                return dict(samples=samples,csq=csq,info_numbers=info_numbers,reference=reference)
            if not line.startswith('#') and line.strip():break
    raise ValueError('This file has no valid VCF #CHROM header. Open a .vcf or .vcf.gz file.')

def resources(build='GRCh38',sfari_path=None):
    genes={};ensembl={}
    with Path(sfari_path or ASSETS/'sfari-genes.csv').open(encoding='utf-8-sig',newline='') as f:
        for r in csv.DictReader(f):
            if 'gene-symbol' not in r:raise ValueError('SFARI CSV needs a gene-symbol column.')
            genes[r['gene-symbol']]=r
            if r.get('ensembl-id'):ensembl[r['ensembl-id'].split('.')[0]]=r
    aliases=dict(l.split() for l in (ASSETS/'to_standard.tsv').read_text().splitlines()) if build=='GRCh38' else {}
    mane=set((ASSETS/'mane_transcripts.txt').read_text().splitlines()) if build=='GRCh38' else set()
    return genes,ensembl,aliases,mane

def annotations(info,alt,alt_index,alts,csq_fields):
    result=[]
    for raw in str(info.get('ANN','')).split(','):
        a=raw.split('|')
        # Match the actual allele. Never copy an annotation across unrelated ALT alleles.
        if len(a)>10 and a[0]==alt:
            result.append(dict(gene=a[3],gene_id=a[4],effect=a[1],impact=a[2],transcript=a[6],hgvs_c=a[9],hgvs_p=a[10],warnings=a[15] if len(a)>15 else ''))
    if result or not csq_fields:return result
    for raw in str(info.get('CSQ','')).split(','):
        a=dict(zip(csq_fields,raw.split('|')))
        num=a.get('ALLELE_NUM')
        if num:
            if num!=str(alt_index):continue
        else:
            # VEP's minimized allele representation: trim shared anchor(s).
            ref=str(info.get('_REF',''));mini=alt;r=ref
            while r and mini and r[0]==mini[0]:r=r[1:];mini=mini[1:]
            while r and mini and r[-1]==mini[-1]:r=r[:-1];mini=mini[:-1]
            if a.get('Allele') not in {alt,mini or '-'}:continue
        result.append(dict(gene=a.get('SYMBOL') or a.get('Gene',''),gene_id=a.get('Gene',''),effect=a.get('Consequence',''),impact=a.get('IMPACT','UNKNOWN'),transcript=a.get('Feature',''),hgvs_c=a.get('HGVSc',''),hgvs_p=a.get('HGVSp',''),warnings=a.get('FLAGS',''),mane_tag=bool(a.get('MANE_SELECT'))))
    return result

SCHEMA='''
CREATE TABLE variants (
 id INTEGER PRIMARY KEY, chromosome TEXT, position INTEGER, ref TEXT, alt TEXT,
 quality REAL, filter TEXT, genotype TEXT, dosage INTEGER, depth REAL, gq REAL,
 allele_depths TEXT, alt_fraction REAL, population_af REAL, assertion TEXT,
 clinvar TEXT, condition TEXT, clinvar_review TEXT, clinvar_allele_id TEXT,
 genes TEXT, effect TEXT, impact TEXT, transcript TEXT, hgvs_c TEXT, hgvs_p TEXT,
 sfari INTEGER, mane INTEGER, review_group TEXT, notes TEXT, annotations TEXT
);
CREATE TABLE variant_genes (
 variant_id INTEGER, gene TEXT, gene_id TEXT, effect TEXT, impact TEXT, transcript TEXT,
 hgvs_c TEXT, hgvs_p TEXT, mane INTEGER, sfari_score TEXT, any_high INTEGER,
 PRIMARY KEY(variant_id,gene)
);
CREATE TABLE metadata (key TEXT PRIMARY KEY,value TEXT);
CREATE TABLE coverage (gene TEXT PRIMARY KEY,bases INTEGER,mean_depth REAL,percent_ge20 REAL,source TEXT);
'''

def allele_value(info,key,index,header,allele_specific=False):
    value=info.get(key)
    if not isinstance(value,str):return value
    number=header['info_numbers'].get(key)
    if number=='A' or allele_specific and ',' in value:
        parts=value.split(',');return parts[index-1] if len(parts)>=index else None
    if number=='R':
        parts=value.split(',');return parts[index] if len(parts)>index else None
    return value

def import_vcf(source,project,sample=None,build='Unknown',progress=None,cancel=None,sfari_path=None):
    source=Path(source).resolve();project=Path(project).resolve();project.mkdir(parents=True,exist_ok=True)
    if build not in ('GRCh38','GRCh37','Unknown'):raise ValueError('Genome build must be GRCh38, GRCh37 or Unknown.')
    header=read_header(source);samples=header['samples']
    if len(samples)>1 and sample is None:raise ValueError('Choose a sample from this multi-sample VCF.')
    sample=sample or (samples[0] if samples else 'Sites only')
    if samples and sample not in samples:raise ValueError(f'Sample {sample!r} is not present in this VCF.')
    sample_col=9+samples.index(sample) if samples else None
    genes,ensembl,aliases,mane_ids=resources(build,sfari_path)
    pending=project/f'.import-{uuid.uuid4().hex}.sqlite';db=sqlite3.connect(pending)
    db.executescript(SCHEMA);db.execute('PRAGMA journal_mode=OFF');db.execute('PRAGMA synchronous=OFF')
    counts=Counter();var_batch=[];gene_batch=[];next_id=0
    try:
        with open_text(source) as f:
            for line_number,line in enumerate(f,1):
                if cancel and cancel.is_set():raise Cancelled('Import cancelled. Source files were left unchanged.')
                if line.startswith('#') or not line.strip():continue
                v=line.rstrip('\r\n').split('\t')
                if len(v)<8:raise ValueError(f'Malformed VCF row at line {line_number}: expected at least 8 columns.')
                try:pos=int(v[1])
                except ValueError:raise ValueError(f'Invalid VCF position at line {line_number}.') from None
                if pos<1:raise ValueError(f'VCF positions must be positive (line {line_number}).')
                if sample_col is not None and len(v)<=sample_col:raise ValueError(f'Missing sample column at line {line_number}.')
                counts['input_records']+=1
                info=split_info(v[7]);info['_REF']=v[3]
                fmt=dict(zip(v[8].split(':'),v[sample_col].split(':'))) if sample_col is not None and len(v)>sample_col else {}
                gt=fmt.get('GT','.')
                alleles=gt.replace('|','/').split('/')
                for alt_index,alt in enumerate(v[4].split(','),1):
                    if alt in ('.','<NON_REF>','*'):
                        counts['placeholder_alleles_skipped']+=1;continue
                    if sample_col is not None and all(a!='.' for a in alleles) and str(alt_index) not in alleles:
                        counts['noncarried_alleles_skipped']+=1;continue
                    dosage=sum(a==str(alt_index) for a in alleles) if all(a!='.' for a in alleles) else None
                    notes=[]
                    if dosage is None:notes.append('Genotype incomplete or unavailable; carrier status unknown')
                    if len(v[4].split(','))>1:notes.append(f'Multiallelic record; this is ALT {alt_index}; original GT retained')
                    chrom=aliases.get(v[0],v[0].removeprefix('chr'))
                    if chrom in ('X','Y'):notes.append('Sex chromosome; ploidy and clinical context require review')
                    anns=annotations(info,alt,alt_index,v[4].split(','),header['csq'])
                    selected={}
                    for a in anns:
                        a['mane']=a['transcript'].split('.')[0] in mane_ids
                        a['mane']=a['mane'] or bool(a.get('mane_tag'))
                        row=genes.get(a['gene']) or ensembl.get(a['gene_id'].split('.')[0])
                        a['sfari_score']=row.get('gene-score','') if row else ''
                        a['sfari']=bool(row)
                        if row:a['gene']=row['gene-symbol']
                        key=a['gene'] or a['gene_id']
                        if not key:continue
                        old=selected.get(key)
                        if old is None or (not a['mane'],IMPACT_ORDER.get(a['impact'],4))<(not old['mane'],IMPACT_ORDER.get(old['impact'],4)):selected[key]=a
                    best=min(selected.values(),key=lambda a:(not a['mane'],IMPACT_ORDER.get(a['impact'],4))) if selected else {}
                    if not anns:counts['unannotated_alleles']+=1;notes.append('No matched ANN/CSQ transcript annotation; gene review unavailable')
                    if any(a['impact']=='HIGH' for a in anns) and best.get('impact')!='HIGH':notes.append('Alternative-transcript HIGH effect; representative annotation differs')
                    dp=numeric(fmt.get('DP'));gq=numeric(fmt.get('GQ'));fraction=None;ad=fmt.get('AD','.')
                    try:
                        depths=[int(x) for x in ad.split(',')]
                        if len(depths)>alt_index and sum(depths):fraction=depths[alt_index]/sum(depths)
                    except ValueError:pass
                    af=numeric(allele_value(info,'POP_AF',alt_index,header,True))
                    # Generic INFO/AF is often a sample-callset statistic, not a population frequency.
                    if af is None:notes.append('Population frequency unavailable; rarity not established')
                    elif not 0<=af<=1:af=None;notes.append('Invalid population frequency ignored')
                    sig=allele_value(info,'CLNSIG',alt_index,header) or ''
                    if len(v[4].split(','))>1 and sig and header['info_numbers'].get('CLNSIG') not in ('A','R'):
                        notes.append('ClinVar annotation may be record-level; verify allele matching')
                    assertion=classify_assertion(sig)
                    good=v[6]=='PASS';support=dp is not None and dp>=10 and gq is not None and gq>=20 and dosage is not None and dosage>0
                    if not support:notes.append('Depth/GQ/genotype support missing or below research thresholds')
                    if fraction is not None and dosage==1 and len(alleles)==2 and not .2<=fraction<=.8:notes.append('Heterozygous allele balance requires review')
                    ambiguous_clinvar=len(v[4].split(','))>1 and bool(sig) and header['info_numbers'].get('CLNSIG') not in ('A','R')
                    if not good:group='Filtered / unfiltered'
                    elif ambiguous_clinvar:group='Support check'
                    elif not support:group='Support check'
                    elif assertion=='Pathogenic / likely pathogenic':group='Clinical assertion'
                    elif assertion=='Conflicting':group='Conflicting assertions'
                    elif assertion=='Benign / likely benign' or af is not None and af>=.01:group='Other'
                    elif best.get('impact')=='HIGH':group='Predicted effect'
                    elif any(a['impact']=='HIGH' for a in anns):group='Transcript check'
                    else:group='Other'
                    next_id+=1
                    row=(next_id,chrom,pos,v[3],alt,numeric(v[5]),v[6],gt,dosage,dp,gq,ad,fraction,af,assertion,str(sig),str(allele_value(info,'CLNDN',alt_index,header) or ''),str(info.get('CLNREVSTAT','')),str(info.get('ALLELEID','')),', '.join(sorted(selected)),best.get('effect','Unannotated'),best.get('impact','UNKNOWN'),best.get('transcript',''),best.get('hgvs_c',''),best.get('hgvs_p',''),int(any(a['sfari'] for a in selected.values())),int(best.get('mane',False)),group,'; '.join(dict.fromkeys(notes)),json.dumps(anns,separators=(',',':')))
                    var_batch.append(row)
                    for gene,a in selected.items():gene_batch.append((next_id,gene,a['gene_id'],a['effect'],a['impact'],a['transcript'],a['hgvs_c'],a['hgvs_p'],int(a['mane']),a['sfari_score'],int(any(x['impact']=='HIGH' and x['gene']==gene for x in anns))))
                    counts['imported_alleles']+=1
                if len(var_batch)>=500:
                    db.executemany('INSERT INTO variants VALUES ('+','.join('?'*30)+')',var_batch);db.executemany('INSERT INTO variant_genes VALUES ('+','.join('?'*11)+')',gene_batch);var_batch.clear();gene_batch.clear();db.commit()
                    if progress:progress(f"Imported {counts['imported_alleles']:,} variant alleles…")
        db.executemany('INSERT INTO variants VALUES ('+','.join('?'*30)+')',var_batch);db.executemany('INSERT INTO variant_genes VALUES ('+','.join('?'*11)+')',gene_batch)
        if not counts['input_records'] and not counts['imported_alleles']:counts['imported_alleles']=0
        db.executescript('CREATE INDEX variant_filter ON variants(filter); CREATE INDEX variant_review ON variants(review_group); CREATE INDEX variant_chrom ON variants(chromosome,position); CREATE INDEX genes_symbol ON variant_genes(gene,variant_id); CREATE INDEX genes_sfari ON variant_genes(sfari_score,variant_id);')
        source_digest=hashlib.sha256()
        with source.open('rb') as source_stream:
            for block in iter(lambda:source_stream.read(1024*1024),b''):source_digest.update(block)
        source_hash=source_digest.hexdigest()
        metadata=dict(schema_version=1,source=str(source),sample=sample,build=build,vcf_reference=header['reference'],imported_at=time.strftime('%Y-%m-%d %H:%M:%S'),counts=dict(counts),sfari_snapshot=json.loads((ASSETS/'resources.json').read_text())['sfari'],representative_annotation='Prefer matching MANE transcript, then predicted impact; all transcript annotations retained',source_sha256=source_hash)
        db.executemany('INSERT INTO metadata VALUES (?,?)',[(k,json.dumps(v)) for k,v in metadata.items()]);db.commit();db.close()
        if cancel and cancel.is_set():raise Cancelled('Import cancelled.')
        os.replace(pending,project/'variants.sqlite')
        (project/'project.json').write_text(json.dumps(metadata,indent=2)+'\n')
        return Project(project)
    except BaseException:
        db.close();pending.unlink(missing_ok=True);raise

class Project:
    def __init__(self,path):
        self.path=Path(path).resolve()
        if not (self.path/'variants.sqlite').exists():raise ValueError('Choose a WES Lab project folder containing variants.sqlite.')
        with self.connect() as db:
            self.meta={r['key']:json.loads(r['value']) for r in db.execute('SELECT * FROM metadata')}
        if self.meta.get('schema_version')!=1:raise ValueError('Unsupported project version.')
    @contextmanager
    def connect(self):
        db=sqlite3.connect(self.path/'variants.sqlite');db.row_factory=sqlite3.Row
        try:
            yield db
            db.commit()
        except BaseException:
            db.rollback();raise
        finally:db.close()
    def where(self,filters=None):
        f=filters or {};where=[];args=[]
        if f.get('pass_only'):where.append("v.filter='PASS'")
        if f.get('sfari_only'):where.append('v.sfari=1')
        if f.get('review_only'):where.append("v.review_group NOT IN ('Other','Filtered / unfiltered')")
        if f.get('group') and f['group']!='All review groups':where.append('v.review_group=?');args.append(f['group'])
        if f.get('impact') and f['impact']!='All effects':where.append('v.impact=?');args.append(f['impact'])
        if f.get('chromosome') and f['chromosome']!='All chromosomes':where.append('v.chromosome=?');args.append(f['chromosome'])
        if f.get('min_dp') is not None:where.append('v.depth>=?');args.append(f['min_dp'])
        if f.get('min_gq') is not None:where.append('v.gq>=?');args.append(f['min_gq'])
        if f.get('query'):
            q='%'+f['query'].strip().upper().replace('!','!!').replace('%','!%').replace('_','!_')+'%'
            where.append("(UPPER(v.genes) LIKE ? ESCAPE '!' OR UPPER(v.hgvs_c) LIKE ? ESCAPE '!' OR UPPER(v.hgvs_p) LIKE ? ESCAPE '!' OR UPPER(v.clinvar) LIKE ? ESCAPE '!' OR UPPER(v.chromosome || ':' || v.position) LIKE ? ESCAPE '!')");args.extend([q]*5)
        return (' WHERE '+' AND '.join(where) if where else ''),args
    def variants(self,filters=None,offset=0,limit=150,sort='position',descending=False):
        fields={'position':"CASE WHEN v.chromosome GLOB '[0-9]*' THEN CAST(v.chromosome AS INTEGER) WHEN v.chromosome='X' THEN 23 WHEN v.chromosome='Y' THEN 24 WHEN v.chromosome IN ('M','MT') THEN 25 ELSE 100 END, v.chromosome, v.position",'gene':'v.genes','depth':'v.depth','gq':'v.gq','impact':"CASE v.impact WHEN 'HIGH' THEN 0 WHEN 'MODERATE' THEN 1 WHEN 'LOW' THEN 2 ELSE 3 END",'review':"CASE v.review_group WHEN 'Clinical assertion' THEN 0 WHEN 'Predicted effect' THEN 1 WHEN 'Conflicting assertions' THEN 2 ELSE 3 END"}
        field=fields.get(sort,fields['position']);where,args=self.where(filters)
        with self.connect() as db:
            count=db.execute('SELECT count(*) FROM variants v'+where,args).fetchone()[0]
            rows=db.execute('SELECT v.* FROM variants v'+where+' ORDER BY '+field+(' DESC' if descending else '')+', v.id LIMIT ? OFFSET ?',args+[limit,offset]).fetchall()
        return count,[dict(r) for r in rows]
    def detail(self,variant_id):
        with self.connect() as db:
            row=db.execute('SELECT * FROM variants WHERE id=?',(variant_id,)).fetchone()
            if row is None:raise ValueError('Variant not found.')
            result=dict(row);result['annotations']=json.loads(result['annotations']);result['gene_annotations']=[dict(x) for x in db.execute('SELECT * FROM variant_genes WHERE variant_id=?',(variant_id,))];return result
    def summary(self):
        with self.connect() as db:
            total=db.execute('SELECT count(*) FROM variants').fetchone()[0]
            passing=db.execute("SELECT count(*) FROM variants WHERE filter='PASS'").fetchone()[0]
            sfari=db.execute('SELECT count(*) FROM variants WHERE sfari=1').fetchone()[0]
            review=db.execute("SELECT count(*) FROM variants WHERE review_group='Clinical assertion'").fetchone()[0]
            groups={r[0]:r[1] for r in db.execute('SELECT review_group,count(*) FROM variants GROUP BY review_group')}
            def distribution(column):return {str(r[0]):r[1] for r in db.execute(f'SELECT {column},count(*) FROM variants GROUP BY {column}')}
            depths={str(r[0]):r[1] for r in db.execute("SELECT CASE WHEN depth IS NULL THEN 'Missing' WHEN depth<10 THEN '<10' WHEN depth<20 THEN '10–19' WHEN depth<50 THEN '20–49' WHEN depth<100 THEN '50–99' ELSE '100+' END,count(*) FROM variants GROUP BY 1")}
            coverage=db.execute('SELECT count(*),sum(percent_ge20>=90) FROM coverage').fetchone()
            types={str(r[0]):r[1] for r in db.execute("SELECT CASE WHEN alt LIKE '<%' OR alt LIKE '%[%' OR alt LIKE '%]%' THEN 'Symbolic / structural' WHEN length(ref)=1 AND length(alt)=1 THEN 'SNV' WHEN length(ref)<length(alt) THEN 'Insertion' WHEN length(ref)>length(alt) THEN 'Deletion' ELSE 'MNV / complex' END,count(*) FROM variants GROUP BY 1")}
            coverage_bins={str(r[0]):r[1] for r in db.execute("SELECT CASE WHEN percent_ge20>=90 THEN '90–100%' WHEN percent_ge20>=50 THEN '50–89%' WHEN percent_ge20>0 THEN 'Below 50%' ELSE 'No bases at 20×' END,count(*) FROM coverage GROUP BY 1")}

            return dict(types=types,coverage_bins=coverage_bins,total=total,passing=passing,sfari=sfari,clinical_assertions=review,groups=groups,impact=distribution('impact'),chromosomes=distribution('chromosome'),assertions=distribution('assertion'),filters=distribution('filter'),depths=depths,coverage_genes=coverage[0],coverage_complete=coverage[1] or 0)
    def import_coverage(self,path,source='Imported table'):
        rows=[]
        with Path(path).open(encoding='utf-8-sig',newline='') as f:
            reader=csv.DictReader(f,delimiter='\t')
            needed={'gene','mean_depth','percent_ge20'}
            if not needed.issubset(reader.fieldnames or []):raise ValueError('Coverage TSV needs gene, mean_depth and percent_ge20 columns.')
            for r in reader:
                mean=numeric(r['mean_depth']);pct=numeric(r['percent_ge20']);bases=numeric(r.get('coding_plus_flank_bases'))
                if not r['gene'] or mean is None or pct is None or mean<0 or not 0<=pct<=100:raise ValueError('Coverage table contains invalid genes, depths or percentages.')
                rows.append((r['gene'],int(bases) if bases is not None else None,mean,pct,source))
        if len({r[0] for r in rows})!=len(rows):raise ValueError('Coverage table has duplicate genes; merge intervals before importing.')
        with self.connect() as db:db.execute('DELETE FROM coverage');db.executemany('INSERT INTO coverage VALUES (?,?,?,?,?)',rows)
        return len(rows)
    def coverage(self,query='',gaps_only=False):
        sql='SELECT * FROM coverage WHERE gene LIKE ?';args=['%'+query.upper()+'%']
        if gaps_only:sql+=' AND percent_ge20<90'
        with self.connect() as db:return [dict(r) for r in db.execute(sql+' ORDER BY percent_ge20,gene',args)]
    def export_csv(self,path,filters=None):
        where,args=self.where(filters)
        with self.connect() as db,Path(path).open('w',encoding='utf-8',newline='') as f:
            cur=db.execute('SELECT * FROM variants v'+where+' ORDER BY chromosome,position,id',args);fields=[c[0] for c in cur.description if c[0]!='annotations'];writer=csv.DictWriter(f,fields);writer.writeheader()
            for r in cur:writer.writerow({k:r[k] for k in fields})
    def export_report(self,path):
        s=self.summary();meta=self.meta;_,rows=self.variants({'group':'Clinical assertion'},limit=100,sort='review')
        def esc(v):return html.escape(str(v if v is not None else '—'))
        cards=''.join(f'<div class="card"><strong>{s[k]:,}</strong><span>{label}</span></div>' for k,label in [('total','Imported variant alleles'),('passing','PASS calls'),('sfari','SFARI overlaps'),('clinical_assertions','Clinical assertions to review')])
        tables=''.join('<tr>'+''.join(f'<td>{esc(v)}</td>' for v in [r['genes'],f"{r['chromosome']}:{r['position']} {r['ref']}>{r['alt']}",r['genotype'],r['depth'],r['gq'],r['hgvs_p'],r['clinvar']])+'</tr>' for r in rows)
        charts=''.join(f'<section><h2>{esc(title)}</h2>{bar_svg(data)}</section>' for title,data in [('Variant types',s['types']),('Predicted effects',s['impact']),('Read depth',s['depths']),('Clinical assertions',s['assertions']),('Measured gene coverage at 20×',s['coverage_bins'])])
        text=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>WES Lab — {esc(meta['sample'])}</title><style>body{{font:15px system-ui;color:#172b32;background:#f5f7f7;max-width:1100px;margin:40px auto;padding:20px}}h1{{font-size:36px}}.cards{{display:flex;gap:14px;flex-wrap:wrap}}.card,section{{background:white;border:1px solid #dce5e4;border-radius:12px;padding:22px;margin:12px 0}}.card strong{{display:block;font-size:28px}}.card span{{color:#547078}}table{{border-collapse:collapse;width:100%;font-size:12px}}td,th{{padding:10px;text-align:left;border-bottom:1px solid #e2e9e8}}svg{{max-width:100%;height:auto}}footer{{color:#567078}}</style><h1>WES Lab / {esc(meta['sample'])}</h1><p>Build: {esc(meta['build'])} · Imported: {esc(meta['imported_at'])}</p><p>Research review. Gene overlap, prediction and ClinVar assertions do not independently establish a diagnosis.</p><div class="cards">{cards}</div>{charts}<section><h2>Clinical assertions — first 100</h2><p>These are database assertions requiring allele, transcript, inheritance and phenotype review. They are not confirmed diagnoses. Quality and support checks use PASS, DP≥10 and GQ≥20; sex-chromosome ploidy remains unresolved.</p><table><tr><th>Gene</th><th>Allele</th><th>GT</th><th>DP</th><th>GQ</th><th>Effect</th><th>ClinVar</th></tr>{tables}</table></section><section><h2>Coverage</h2><p>{s['coverage_complete']:,} of {s['coverage_genes']:,} measured genes have at least 90% of supplied coding/flank bases at 20×. Coverage results come from the attached table; coverage cannot be inferred from a variant-only VCF.</p></section><footer>Source: {esc(meta['source'])}<br>SHA256: {esc(meta['source_sha256'])}<br>SFARI snapshot: {esc(meta['sfari_snapshot']['local_snapshot'])}. Population-frequency absence is not evidence of rarity. No alignment, CNV, repeat expansion or clinical confirmation is implied by VCF import.</footer></html>'''
        Path(path).write_text(text,encoding='utf-8')

def bar_svg(data,width=760):
    items=list(data.items());height=max(100,45*len(items)+20);maxv=max(data.values(),default=1) or 1
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="Variant count bar chart">']
    for i,(label,value) in enumerate(items):
        y=20+i*45;w=(width-330)*value/maxv
        parts.append(f'<text x="0" y="{y+15}" fill="#46616a" font-family="sans-serif" font-size="13">{html.escape(str(label))}</text><rect x="220" y="{y}" width="{w:.2f}" height="22" rx="4" fill="#158f86"/><text x="{230+w:.2f}" y="{y+15}" fill="#172b32" font-family="sans-serif" font-size="12">{value:,}</text>')
    return ''.join(parts)+'</svg>'
