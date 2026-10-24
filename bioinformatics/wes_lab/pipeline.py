"""Sequential local calling/annotation for an analysis-ready, compatible GRCh38 BAM.

Does not align FASTQs, mark duplicates, recalibrate reads or make diagnoses.
The installed accession-based GRCh38 resource bundle is required.
"""
from __future__ import annotations
import collections,contextlib,hashlib,json,os,re,signal,subprocess,threading,time
from pathlib import Path
from .data import Cancelled,open_text,read_header
ROOT=Path(__file__).resolve().parents[1]

def bundle():
    return {'samtools':ROOT/'tools/samtools-env/bin/samtools','bcftools':ROOT/'tools/variants-env/bin/bcftools','gatk':ROOT/'tools/gatk-4.7.0.0/gatk','snpeff':ROOT/'tools/variants-env/bin/snpEff','reference':ROOT/'data/reference/GRCh38.fna','targets':ROOT/'data/reference/wes.exons_padded100.bed','clinvar':ROOT/'data/reference/clinvar.GRCh38.matched.vcf.gz','frequency':ROOT/'data/reference/1000G.af_available.vcf.gz','standard':ROOT/'data/reference/to_standard.tsv','accessions':ROOT/'data/reference/to_accessions.tsv','annotation_bed':ROOT/'data/reference/ensembl115.bed','snpeff_data':ROOT/'data/snpeff','predictor':ROOT/'data/snpeff/GRCh38.115/snpEffectPredictor.bin'}

def probe_bundle():
    return {k:str(p) for k,p in bundle().items() if not p.exists()}

class Pipeline:
    def __init__(self,source,output,mode='bam',progress=None,cancel=None):
        self.source=Path(source).resolve();self.out=Path(output).resolve();self.mode=mode
        if mode not in ('bam','vcf'):raise ValueError('Unsupported workflow.')
        self.progress=progress or (lambda _:None);self.cancel=cancel or threading.Event();self.process=None;self.paths=bundle()
        self.out.mkdir(parents=True,exist_ok=True);(self.out/'parts').mkdir(exist_ok=True)
    def status(self,message):
        self.progress(message);(self.out/'STATUS.txt').write_text(time.strftime('%Y-%m-%d %H:%M:%S')+' '+message+'\n')
    def execute(self,args,stdout=None):
        if self.cancel.is_set():raise Cancelled('Workflow stopped; completed parts are saved.')
        args=[str(x) for x in args]
        with (self.out/'workflow.log').open('ab') as log:
            log.write(('\nCOMMAND '+repr(args)+'\n').encode());log.flush()
            with contextlib.ExitStack() as stack:
                dest=stack.enter_context(Path(stdout).open('wb')) if stdout else log
                self.process=subprocess.Popen(args,stdout=dest,stderr=log,start_new_session=True)
                try:
                    while self.process.poll() is None:
                        if self.cancel.wait(.3):
                            with contextlib.suppress(ProcessLookupError):os.killpg(self.process.pid,signal.SIGTERM)
                            try:self.process.wait(timeout=5)
                            except subprocess.TimeoutExpired:
                                with contextlib.suppress(ProcessLookupError):os.killpg(self.process.pid,signal.SIGKILL)
                                self.process.wait()
                            raise Cancelled('Workflow stopped; completed parts are saved.')
                    if self.process.returncode:raise RuntimeError(f'{Path(args[0]).name} failed. See {self.out / "workflow.log"}.')
                finally:self.process=None
    def bcf(self,*args,stdout=None):self.execute([self.paths['bcftools'],*args],stdout)
    def index(self,path):self.bcf('index','-f','-t',path)
    def publish(self,partial,target):
        self.index(partial);os.replace(partial,target);os.replace(str(partial)+'.tbi',str(target)+'.tbi')
    def preflight(self):
        missing=probe_bundle()
        if missing:raise ValueError('Missing local resources: '+', '.join(missing)+'. See WES_LAB.md for setup requirements.')
        if not self.source.is_file():raise ValueError('Input file is missing.')
        if not Path(str(self.paths['reference'])+'.fai').exists() or not self.paths['reference'].with_suffix('.dict').exists():raise ValueError('Reference FASTA index or dictionary is missing.')
        names={l.split('\t')[0] for l in Path(str(self.paths['reference'])+'.fai').read_text().splitlines()}
        if self.mode=='bam':
            if not any(p.exists() for p in (Path(str(self.source)+'.bai'),self.source.with_suffix('.bai'),Path(str(self.source)+'.csi'))):raise ValueError('BAM index is missing. Supply a sorted, indexed, analysis-ready BAM.')
            result=subprocess.run([str(self.paths['samtools']),'view','-H',str(self.source)],capture_output=True,text=True,check=True)
            header=result.stdout
            if not re.search(r'^@HD\t.*SO:coordinate',header,re.M):raise ValueError('BAM header must declare coordinate sorting.')
            sq={dict(x.split(':',1) for x in l.split('\t')[1:] if ':' in x).get('SN'):int(dict(x.split(':',1) for x in l.split('\t')[1:] if ':' in x).get('LN','0')) for l in header.splitlines() if l.startswith('@SQ\t')}
            reference_lengths={p[0]:int(p[1]) for l in Path(str(self.paths['reference'])+'.fai').read_text().splitlines() if (p:=l.split('\t'))}
            targeted={l.split('\t')[0] for l in self.paths['targets'].read_text().splitlines()}
            if not targeted.issubset(sq) or any(sq[c]!=reference_lengths[c] for c in targeted):raise ValueError('BAM contig names/lengths do not match the installed GRCh38 reference and targets. Use a compatible BAM; do not rename contigs across genome builds.')
            samples={x[3:] for l in header.splitlines() if l.startswith('@RG\t') for x in l.split('\t') if x.startswith('SM:')}
            if len(samples)!=1:raise ValueError('BAM must contain read groups for exactly one named sample.')
            self.sample=next(iter(samples))
            self.execute([self.paths['samtools'],'quickcheck','-v',self.source])
        else:
            header=read_header(self.source)
            if len(header['samples'])!=1:raise ValueError('Annotation workflow expects a single-sample VCF. Import multi-sample VCFs directly into the viewer instead.')
            self.sample=header['samples'][0]
            source_contigs=set()
            with open_text(self.source) as f:
                for line in f:
                    if not line.startswith('#') and line.strip():source_contigs.add(line.split('\t',1)[0])
            primary=set(dict(l.split() for l in self.paths['standard'].read_text().splitlines()))
            if not source_contigs.issubset(primary):raise ValueError('Annotation workflow expects nuclear primary accession-named GRCh38 variants matching the installed bundle. The viewer can open other builds; this workflow cannot annotate them.')
        # Fingerprint configuration and resource stats, avoiding stale checkpoints.
        config={'version':1,'implementation_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'mode':self.mode,'source':str(self.source),'source_size':self.source.stat().st_size,'source_mtime_ns':self.source.stat().st_mtime_ns,'sample':self.sample,'resources':{k:[str(p),p.stat().st_size,p.stat().st_mtime_ns] for k,p in self.paths.items()},'ploidy':2,'annotation_window':10_000_000,'interval_padding':100,'java_annotation':'3g_serial_1cpu'}
        manifest=self.out/'configuration.json'
        if manifest.exists() and json.loads(manifest.read_text())!=config:raise ValueError('Input or resources changed. Choose a new workflow folder to preserve completed outputs.')
        manifest.write_text(json.dumps(config,indent=2)+'\n')
    def run(self):
        # Linux process lock, released even if the process crashes.
        import fcntl
        with (self.out/'workflow.lock').open('w') as lock:
            try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError:raise ValueError('This workflow is already running.') from None
            try:
                self.status('Checking inputs and the local GRCh38 resource bundle');self.preflight()
                source=self.call_bam() if self.mode=='bam' else self.source
                result=self.annotate(source)
                self.status('COMPLETE: annotated variants ready to import');return result
            except Cancelled:
                self.status('STOPPED: restart with the same input and folder to resume');raise
            except Exception:
                self.status('FAILED: see workflow.log; completed parts retained');raise
    def call_bam(self):
        target_groups=collections.defaultdict(list)
        for line in self.paths['targets'].read_text().splitlines():target_groups[line.split('\t')[0]].append(line)
        completed=[]
        for i,(chrom,lines) in enumerate(target_groups.items(),1):
            if self.cancel.is_set():raise Cancelled()
            base=self.out/'parts'/f'call_{chrom}';called=Path(str(base)+'.vcf.gz')
            if not Path(str(base)+'.complete').exists():
                self.status(f'Calling chromosome {i}/{len(target_groups)}: {chrom} (diploid; X/Y provisional)')
                bed=Path(str(base)+'.bed');bed.write_text('\n'.join(lines)+'\n');raw=Path(str(base)+'.raw.vcf.gz')
                self.execute([self.paths['gatk'],'--java-options','-Xmx2g -XX:ParallelGCThreads=1', 'HaplotypeCaller','-R',self.paths['reference'],'-I',self.source,'-L',bed,'--interval-padding','100','--sample-ploidy','2','--native-pair-hmm-threads','1','-O',raw])
                self.bcf('view','-R',bed,'--regions-overlap','0',raw,'-Oz','-o',called);self.index(called);Path(str(base)+'.complete').touch()
            completed.append(called)
        merged=self.out/'called.vcf.gz';filtered=self.out/'filtered.vcf.gz'
        if not (self.out/'filtering.complete').exists():
            self.status('Merging calls, normalizing and applying exploratory quality filters')
            listing=self.out/'calling.list';listing.write_text(''.join(str(p)+'\n' for p in completed));self.bcf('concat','-f',listing,'-Oz','-o',merged)
            normal=self.out/'normalized.vcf.gz';unsorted=self.out/'normalized.unsorted.vcf.gz'
            self.bcf('norm','-f',self.paths['reference'],'-c','e','-m','-any',merged,'-Oz','-o',unsorted)
            self.bcf('sort','-m','128M','-T',self.out/'sort_tmp','-Oz','-o',normal,unsorted)
            snps=self.out/'snp_filtered.vcf.gz';indels=self.out/'indel_filtered.vcf.gz'
            self.bcf('filter','-m','+','-s','SNPQuality','-e','TYPE="snp" && (QUAL<30 || INFO/QD<2 || INFO/FS>60 || INFO/SOR>3 || INFO/MQ<40 || INFO/MQRankSum< -12.5 || INFO/ReadPosRankSum< -8)',normal,'-Oz','-o',snps)
            self.bcf('filter','-m','+','-s','IndelQuality','-e','TYPE="indel" && (QUAL<30 || INFO/QD<2 || INFO/FS>200 || INFO/ReadPosRankSum< -20)',snps,'-Oz','-o',indels)
            self.bcf('filter','-m','+','-s','LowGenotypeSupport','-e','FMT/DP="." || FMT/GQ="." || FMT/DP<10 || FMT/GQ<20',indels,'-Oz','-o',filtered);self.index(filtered);(self.out/'filtering.complete').touch()
        return filtered
    def annotate(self,source):
        # Normalize even imported VCFs to biallelic, reference-checked, indexed input.
        prepared=self.out/'annotation_input.vcf.gz'
        if not (self.out/'annotation_input.complete').exists():
            unsorted=self.out/'annotation_input.unsorted.vcf.gz'
            self.bcf('norm','-f',self.paths['reference'],'-c','e','-m','-any',source,'-Oz','-o',unsorted)
            self.bcf('sort','-m','128M','-T',self.out/'sort_tmp','-Oz','-o',prepared,unsorted);self.index(prepared);(self.out/'annotation_input.complete').touch()
        windows=collections.OrderedDict();input_keys=collections.Counter()
        with open_text(prepared) as f:
            for line in f:
                if line.startswith('#'):continue
                v=line.split('\t');key=(v[0],v[1],v[3],v[4]);input_keys[key]+=1
                start=((int(v[1])-1)//10_000_000)*10_000_000+1;windows.setdefault((v[0],start),0);windows[(v[0],start)]+=1
        if not windows:raise ValueError('No variants were called in the chosen intervals.')
        aliases=dict(l.split() for l in self.paths['standard'].read_text().splitlines());transcripts=collections.defaultdict(set)
        with self.paths['annotation_bed'].open() as f:
            for line in f:
                c,a,b,feature=line.rstrip().split('\t')[:4]
                if feature.startswith('Transcript_'):transcripts[c].add(feature[11:])
        completed=[]
        for i,((chrom,start),expected) in enumerate(windows.items(),1):
            if self.cancel.is_set():raise Cancelled()
            if chrom not in aliases or not transcripts[aliases[chrom]]:raise ValueError(f'No transcript/chromosome mapping for {chrom}.')
            base=self.out/'parts'/f'ann_{chrom}_{start}';final=Path(str(base)+'.vcf.gz');marker=Path(str(base)+'.complete')
            if not marker.exists():
                self.status(f'Annotating section {i}/{len(windows)}: {aliases[chrom]}:{start:,}–{start+9_999_999:,}')
                ids=Path(str(base)+'.transcripts.txt');ids.write_text('\n'.join(sorted(transcripts[aliases[chrom]]))+'\n')
                raw=Path(str(base)+'.input.vcf.gz');clin=Path(str(base)+'.clinvar.vcf.gz');pop=Path(str(base)+'.population.vcf.gz');standard=Path(str(base)+'.standard.vcf');effects=Path(str(base)+'.effects.vcf');partial=Path(str(base)+'.partial.vcf.gz')
                self.bcf('view','-r',f'{chrom}:{start}-{start+9_999_999}','--regions-overlap','0',prepared,'-Oz','-o',raw);self.index(raw)
                self.bcf('annotate','-a',self.paths['clinvar'],'-c','INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT,INFO/ALLELEID',raw,'-Oz','-o',clin);self.index(clin)
                self.bcf('annotate','-a',self.paths['frequency'],'-c','INFO/POP_AF:=INFO/AF',clin,'-Oz','-o',pop)
                self.bcf('annotate','--rename-chrs',self.paths['standard'],pop,'-Ov','-o',standard)
                self.execute([self.paths['snpeff'],'-Xmx3g','-XX:+UseSerialGC','-XX:ActiveProcessorCount=1','-onlyTr',ids,'-dataDir',self.paths['snpeff_data'],'-noStats','-noLog','-nodownload','GRCh38.115',standard],stdout=effects)
                self.bcf('annotate','--rename-chrs',self.paths['accessions'],effects,'-Oz','-o',partial);self.publish(partial,final)
                observed=sum(1 for l in open_text(final) if not l.startswith('#'))
                if observed!=expected:raise ValueError('Annotation record count changed; section not marked complete.')
                marker.write_text(str(expected)+'\n')
            completed.append(final)
        self.status('Merging sections and checking every chromosome/position/REF/ALT')
        listing=self.out/'annotation.list';listing.write_text(''.join(str(p)+'\n' for p in completed));partial=self.out/'annotated.partial.vcf.gz';final=self.out/'annotated.vcf.gz'
        self.bcf('concat','-f',listing,'-Oz','-o',partial)
        observed=collections.Counter()
        with open_text(partial) as f:
            for l in f:
                if not l.startswith('#'):
                    v=l.split('\t');observed[(v[0],v[1],v[3],v[4])]+=1
        if input_keys!=observed:raise ValueError('Merged variants differ from the normalized input; final file not published.')
        self.publish(partial,final);return final
