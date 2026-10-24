"""Local IGV desktop controller for the reviewed 571566 BAM/VCF bundle."""
from __future__ import annotations
import argparse,csv,json,queue,socket,subprocess,threading,time
from pathlib import Path
import xml.etree.ElementTree as ET
import tkinter as tk
from tkinter import ttk,messagebox
from .people import atomic_json,lock_file

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/wgs/igv_review'
QUALITY=ROOT/'results/wgs/quality_review'
LABELS={('17',74863914):'P/LP — FDXR',('2',218889989):'P/LP — WNT10A',('12',6034812):'Pathogenic, recessive — VWF',('14',60979582):'Conflicting — TRMT5',('9',128945903):'Unclassified — DOLK',('1',225240648):'Single research assertion — DNAH14',('19',1037757):'Unreliable call — CNN2',('19',1037767):'Unreliable call — CNN2'}

def locus(chrom,position,flank=100):
    aliases=dict(l.split() for l in (ROOT/'data/reference/to_standard.tsv').read_text().splitlines())
    reverse={v:k for k,v in aliases.items()}
    chrom=str(chrom).removeprefix('chr');chrom=reverse.get(chrom,chrom)
    if chrom not in aliases:raise ValueError('Chromosome is absent from the configured GRCh38 reference.')
    position=int(position)
    if position<1:raise ValueError('Position must be positive.')
    return f'{chrom}:{max(1,position-flank)}-{position+flank}'

def prepare_session(chrom='17',position=74863914):
    OUT.mkdir(parents=True,exist_ok=True)
    required=[ROOT/'data/reference/GRCh38.fna',ROOT/'data/reference/GRCh38.fna.fai',ROOT/'results/wgs/571566.recal.bam',ROOT/'results/wgs/571566.recal.bai',QUALITY/'571566.quality_accepted.vcf.gz',QUALITY/'571566.quality_accepted.vcf.gz.tbi',QUALITY/'571566.needs_quality_review.vcf.gz',QUALITY/'571566.needs_quality_review.vcf.gz.tbi']
    for path in required:
        if not path.is_file():raise ValueError(f'Required local file is missing: {path}')
    aliases=dict(l.split() for l in (ROOT/'data/reference/to_standard.tsv').read_text().splitlines());reverse={v:k for k,v in aliases.items()}
    with (QUALITY/'report_findings_quality.tsv').open() as f,(OUT/'review_sites.bed').open('w') as bed:
        for r in csv.DictReader(f,delimiter='\t'):
            pos=int(r['position']);bed.write(f"{reverse[r['chromosome']]}\t{pos-1}\t{pos-1+len(r['ref'])}\t{r['genes']}:{r['quality_status']}\n")
    root=ET.Element('Session',genome=str(required[0]),locus=locus(chrom,position),version='8')
    resources=ET.SubElement(root,'Resources')
    for path,name in [(required[2],'571566 — recalibrated reads'),(required[4],'Quality screen — retained'),(required[6],'Quality screen — needs review'),(OUT/'review_sites.bed','Report findings')]:
        ET.SubElement(resources,'Resource',path=str(path),name=name)
    path=OUT/'571566.xml';ET.ElementTree(root).write(path,encoding='utf-8',xml_declaration=True)
    return path

def send_command(port,command):
    if '\n' in command or '\r' in command:raise ValueError('Invalid IGV command')
    with socket.create_connection(('127.0.0.1',port),timeout=3) as conn:
        conn.settimeout(15);conn.sendall((command+'\n').encode())
        result=conn.recv(65536).decode().strip()
        if not result or result.lower().startswith(('error','unknown')):raise ValueError('IGV: '+result)
        return result

class Review(tk.Tk):
    def __init__(self,chrom='17',position=74863914):
        super().__init__();self.title('WES Lab — IGV variant review / 571566');self.geometry('1250x760')
        self.process=None;self.port=None;self.log=None;self.events=queue.Queue();self.current=None;self.pending=None
        self.session=prepare_session(chrom,position)
        with (QUALITY/'variant_quality.tsv').open() as f:self.rows=list(csv.DictReader(f,delimiter='\t'))
        self.notes_path=OUT/'visual_review.json'
        ttk.Label(self,text='Inspect the reads behind each variant',font=('DejaVu Sans',18,'bold')).pack(anchor='w',padx=20,pady=(20,5))
        ttk.Label(self,text='Sample 571566 · local GRCh38 reference · BAM alignments + retained/review VCF tracks · IGV memory limit 2 GB').pack(anchor='w',padx=20)
        controls=ttk.Frame(self);controls.pack(fill='x',padx=20,pady=15)
        self.search=tk.StringVar();ttk.Entry(controls,textvariable=self.search,width=25).pack(side='left');self.search.trace_add('write',lambda *_:self.refresh())
        self.scope=tk.StringVar(value='Report findings');cb=ttk.Combobox(controls,textvariable=self.scope,values=['Report findings','All variants','Quality retained','Needs quality review'],state='readonly',width=22);cb.pack(side='left',padx=10);cb.bind('<<ComboboxSelected>>',lambda _:self.refresh())
        ttk.Button(controls,text='Open selected in IGV',command=self.open_selected).pack(side='right')
        cols=['gene','locus','gt','dp','gq','quality','classification','review']
        frame=ttk.Frame(self);frame.pack(fill='both',expand=True,padx=20)
        self.tree=ttk.Treeview(frame,columns=cols,show='headings',selectmode='browse');scroll=ttk.Scrollbar(frame,orient='vertical',command=self.tree.yview);self.tree.configure(yscrollcommand=scroll.set)
        for c,w in zip(cols,[160,140,60,50,50,120,230,150]):self.tree.heading(c,text=c.title());self.tree.column(c,width=w,minwidth=40)
        self.tree.pack(side='left',fill='both',expand=True);scroll.pack(side='right',fill='y');self.tree.bind('<<TreeviewSelect>>',self.selected);self.tree.bind('<Double-1>',lambda _:self.open_selected())
        self.detail=tk.StringVar();ttk.Label(self,textvariable=self.detail,wraplength=1150).pack(fill='x',padx=20,pady=10)
        edit=ttk.Frame(self);edit.pack(fill='x',padx=20,pady=8)
        self.verdict=tk.StringVar(value='Unreviewed');ttk.Combobox(edit,textvariable=self.verdict,values=['Unreviewed','Reads support call','Suspect artifact','Needs confirmation'],state='readonly',width=22).pack(side='left')
        self.note=tk.StringVar();ttk.Entry(edit,textvariable=self.note,width=65).pack(side='left',padx=10);ttk.Button(edit,text='Save visual review',command=self.save).pack(side='left')
        ttk.Button(edit,text='Export reviews TSV',command=self.export).pack(side='right')
        self.status=tk.StringVar();ttk.Label(self,textvariable=self.status).pack(fill='x',padx=20,pady=10)
        ttk.Label(self,text='Read-level review is not a clinical diagnosis. Database classifications are retained review labels; saving a visual review does not change pathogenicity.').pack(anchor='w',padx=20,pady=(0,15))
        self.refresh();self.after(100,self.poll)
        match=next((i for i,r in self.visible.items() if r['chromosome']==chrom and int(r['position'])==int(position)),None)
        if match:self.tree.selection_set(match);self.tree.see(match)
        self.after(500,self.open_selected)
    def notes(self):
        if not self.notes_path.exists():return {}
        return json.loads(self.notes_path.read_text())
    @staticmethod
    def key(r):return ':'.join(r[k] for k in ('chromosome','position','ref','alt'))
    def refresh(self):
        query=self.search.get().lower().strip();scope=self.scope.get();notes=self.notes();matched=[]
        for i,r in enumerate(self.rows):
            if scope=='Report findings' and (r['chromosome'],int(r['position'])) not in LABELS:continue
            if scope=='Quality retained' and r['quality_status']!='ACCEPTED':continue
            if scope=='Needs quality review' and r['quality_status']!='NEEDS_REVIEW':continue
            if query and query not in (' '.join([r['genes'],self.key(r),LABELS.get((r['chromosome'],int(r['position'])), '')])).lower():continue
            matched.append((str(i),r))
        self.tree.delete(*self.tree.get_children());self.visible=dict(matched[:500])
        for i,r in self.visible.items():self.tree.insert('','end',iid=i,values=(r['genes'],f"{r['chromosome']}:{r['position']}",r['GT'],r['DP'],r['GQ'],r['quality_status'],LABELS.get((r['chromosome'],int(r['position'])),'Not curated'),notes.get(self.key(r),{}).get('status','Unreviewed')))
        self.status.set(f'{len(matched):,} matches · showing first {min(500,len(matched)):,}; narrow by gene or chromosome:position.')
    def selected(self,_=None):
        ids=self.tree.selection()
        self.current=self.visible.get(ids[0]) if ids else None
        if self.current:
            r=self.current;n=self.notes().get(self.key(r),{});self.verdict.set(n.get('status','Unreviewed'));self.note.set(n.get('note',''))
            self.detail.set(f"{self.key(r)} · AD {r['AD']} · ALT fraction {r['ALT_fraction']} · MQ {r['MQ']} · Quality reasons: {r['failure_reasons'] or 'none'} · Warnings: {r['warnings'] or 'none'}")
    def save(self):
        if not self.current:return
        try:
            with lock_file(OUT/'.reviews.lock'):
                notes=self.notes();notes[self.key(self.current)]={'status':self.verdict.get(),'note':self.note.get(),'updated_at':time.strftime('%Y-%m-%d %H:%M:%S')};atomic_json(self.notes_path,notes)
            self.refresh();self.status.set('Visual review saved locally.')
        except (OSError,ValueError) as e:messagebox.showerror('Could not save',str(e))
    def export(self):
        with (OUT/'visual_review.tsv').open('w') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(['allele','review_status','note','updated_at'])
            for key,n in self.notes().items():w.writerow([key,n['status'],n['note'],n['updated_at']])
        self.status.set(f'Exported {OUT / "visual_review.tsv"}')
    def open_selected(self):
        self.selected()
        if not self.current:return
        if self.pending:return
        target=locus(self.current['chromosome'],self.current['position'])
        self.pending=True;self.status.set('Opening local reads in IGV…')
        threading.Thread(target=self.navigate,args=(target,),daemon=True).start()
    def navigate(self,target):
        try:
            if self.process is None or self.process.poll() is not None:
                with socket.socket() as s:s.bind(('127.0.0.1',0));self.port=s.getsockname()[1]
                home=OUT/'profile';home.mkdir(exist_ok=True);igv=ROOT/'tools/IGV_2.19.8'
                self.log=(OUT/'igv.log').open('w')
                java=ROOT/'tools/IGV_Linux_2.19.8/jdk-21/bin/java'
                if not java.exists():raise ValueError('Install the bundled graphical Java runtime before opening IGV.')
                self.process=subprocess.Popen([str(java),'-Xmx2g',f'-Duser.home={home}','--module-path='+str(igv/'lib'),'@'+str(igv/'igv.args'),'--module=org.igv/org.broad.igv.ui.Main','-p',str(self.port),str(self.session)],cwd=ROOT,stdout=self.log,stderr=subprocess.STDOUT)
                # IGV opens its command listener before the initial session has
                # loaded a genome. A premature goto can trigger an AWT error.
                deadline=time.monotonic()+90
                while True:
                    if self.process.poll() is not None:raise ValueError('IGV exited; see igv.log')
                    log=(OUT/'igv.log').read_text(errors='replace')
                    if '[TrackLoader] Loading resource:  '+str(OUT/'review_sites.bed') in log:
                        time.sleep(2);break
                    if time.monotonic()>deadline:raise ValueError('IGV session did not finish loading; see igv.log')
                    time.sleep(.5)
            deadline=time.monotonic()+90
            while True:
                if self.process.poll() is not None:raise ValueError(f'IGV exited; see {OUT / "igv.log"}')
                try:send_command(self.port,'goto '+target);break
                except (OSError,ValueError):
                    if time.monotonic()>deadline:raise ValueError('IGV did not become ready; see igv.log')
                    time.sleep(.5)
            self.events.put(('ok','IGV opened at '+target+' · inspect read bases, strands, alignment and coverage.'))
        except Exception as e:self.events.put(('error',str(e)))
    def poll(self):
        try:
            while True:
                kind,text=self.events.get_nowait();self.pending=False;self.status.set(text)
                if kind=='error':messagebox.showerror('IGV',text)
        except queue.Empty:pass
        self.after(100,self.poll)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--chrom',default='17');parser.add_argument('--position',type=int,default=74863914);args=parser.parse_args()
    Review(args.chrom,args.position).mainloop()

if __name__=='__main__':main()
