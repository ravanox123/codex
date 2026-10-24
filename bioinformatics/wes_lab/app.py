"""Native Tk desktop application. Long-running work stays off the GUI thread."""
from __future__ import annotations
import argparse,csv,json,math,queue,threading,time,traceback,webbrowser,uuid
from pathlib import Path
import tkinter as tk
from tkinter import ttk,filedialog,messagebox
from tkinter.scrolledtext import ScrolledText
from urllib.parse import quote
from .data import Project,import_vcf,read_header,bar_svg,Cancelled
from .pipeline import Pipeline,probe_bundle,ROOT
from .people import PersonProject,find_analysis_owner
from .person_ui import PeopleUI

BG='#f3f6f5';PAPER='#ffffff';INK='#18373e';MUTED='#637c80';TEAL='#158f83';NAV='#122d34';BORDER='#dce6e3';AMBER='#b97829'
COLORS=['#158f83','#539eae','#e2b86d','#a18abd','#799c86','#d88f7c']

class Chart(tk.Frame):
    def __init__(self,parent,title,subtitle='',data=None,on_select=None):
        super().__init__(parent,bg=PAPER,highlightbackground=BORDER,highlightthickness=1)
        self.data=data or {};self.on_select=on_select
        heading=tk.Frame(self,bg=PAPER);heading.pack(fill='x',padx=20,pady=(17,0))
        tk.Label(heading,text=title,bg=PAPER,fg=INK,font=('DejaVu Sans',12,'bold')).pack(side='left')
        tk.Button(heading,text='Save graph',command=self.export,bg=PAPER,fg=TEAL,relief='flat',cursor='hand2',font=('DejaVu Sans',9)).pack(side='right')
        tk.Label(self,text=subtitle,bg=PAPER,fg=MUTED,font=('DejaVu Sans',9),anchor='w').pack(fill='x',padx=20,pady=(5,0))
        self.canvas=tk.Canvas(self,bg=PAPER,highlightthickness=0,height=max(200,len(self.data)*29+30))
        self.canvas.pack(fill='both',expand=True,padx=16,pady=12);self.canvas.bind('<Configure>',lambda _:self.draw())
        self.tip=tk.Label(self,text='',bg=PAPER,fg=MUTED,font=('DejaVu Sans',9),anchor='w');self.tip.pack(fill='x',padx=20,pady=(0,12))
    def draw(self):
        c=self.canvas;c.delete('all');w=c.winfo_width();h=c.winfo_height()
        if not self.data or not any(self.data.values()):c.create_text(w/2,h/2,text='No data to display',fill=MUTED,font=('DejaVu Sans',11));return
        items=list(self.data.items());maxv=max(v for _,v in items) or 1;left=min(180,w*.38);available=max(20,w-left-70);step=max(21,min(34,(h-12)/max(1,len(items))))
        for i,(label,value) in enumerate(items):
            y=8+i*step;bar=available*value/maxv;color=COLORS[i%len(COLORS)]
            c.create_text(3,y+9,text=str(label),anchor='w',fill=MUTED,font=('DejaVu Sans',9),width=left-8)
            c.create_rectangle(left,y,left+available,y+17,fill='#edf3f1',outline='')
            item=c.create_rectangle(left,y,left+max(bar,1),y+17,fill=color,outline='')
            c.create_text(left+bar+7,y+9,text=f'{value:,}',anchor='w',fill=INK,font=('DejaVu Sans',9))
            c.tag_bind(item,'<Enter>',lambda _,label=label,value=value:self.tip.configure(text=f'{label}: {value:,}'+(' · click to explore' if self.on_select else '')))
            if self.on_select:c.tag_bind(item,'<Button-1>',lambda _,label=label:self.on_select(label))
    def export(self):
        path=filedialog.asksaveasfilename(title='Save graph',defaultextension='.svg',filetypes=[('SVG image','*.svg')])
        if path:
            try:Path(path).write_text(bar_svg(self.data),encoding='utf-8')
            except OSError as e:messagebox.showerror('Could not save graph',str(e))

class ScrollPage(tk.Frame):
    def __init__(self,parent):
        super().__init__(parent,bg=BG)
        canvas=tk.Canvas(self,bg=BG,highlightthickness=0);bar=ttk.Scrollbar(self,orient='vertical',command=canvas.yview);canvas.configure(yscrollcommand=bar.set)
        bar.pack(side='right',fill='y');canvas.pack(side='left',fill='both',expand=True)
        self.body=tk.Frame(canvas,bg=BG);window=canvas.create_window(0,0,window=self.body,anchor='nw')
        self.body.bind('<Configure>',lambda _:canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.bind('<Configure>',lambda e:canvas.itemconfigure(window,width=e.width))
        def wheel(e):
            if not self.winfo_exists() or not str(e.widget).startswith(str(self)):return
            if e.num==4:canvas.yview_scroll(-3,'units')
            elif e.num==5:canvas.yview_scroll(3,'units')
            elif e.delta:canvas.yview_scroll(int(-e.delta/120)*3,'units')
        self.bind_all('<Button-4>',wheel,add='+');self.bind_all('<Button-5>',wheel,add='+');self.bind_all('<MouseWheel>',wheel,add='+')

class App(PeopleUI,tk.Tk):
    def __init__(self,example=False,project=None,restore=True):
        super().__init__();self.title('WES Lab — local exome research');self.geometry('1320x850');self.minsize(1050,700);self.configure(bg=BG)
        self.person=None;self.project=None;self.page='Projects';self.busy=False;self.cancel=threading.Event();self.events=queue.Queue();self.worker=None;self.offset=0;self.page_size=150;self.debounce=None;self.summary_cache={}
        self.workspace=ROOT/'results/wes_lab';self.workspace.mkdir(parents=True,exist_ok=True);self.settings_path=self.workspace/'settings.json';self.people_root=self.workspace/'people';self.people_root.mkdir(exist_ok=True)
        try:self.settings=json.loads(self.settings_path.read_text())
        except (OSError,ValueError):self.settings={'recent':[]}
        self.configure_styles();self.make_shell();self.show_page('Projects');self.after(100,self.poll);self.protocol('WM_DELETE_WINDOW',self.close)
        if project:self.after(200,lambda:self.open_project(Path(project)))
        elif example:self.after(200,self.open_example)
        elif restore and self.settings.get('last_person') and (Path(self.settings['last_person'])/'person.json').exists():self.after(200,lambda:self.open_person(self.settings['last_person']))
        elif restore and self.settings.get('recent') and (Path(self.settings['recent'][0])/'variants.sqlite').exists():self.after(200,lambda:self.open_project(self.settings['recent'][0]))
    def configure_styles(self):
        s=ttk.Style(self);s.theme_use('clam');s.configure('.',font=('DejaVu Sans',10),background=BG,foreground=INK)
        s.configure('TButton',padding=(13,8),borderwidth=0,background='#e4eeea',foreground=INK);s.map('TButton',background=[('active','#d6e6e0'),('disabled','#edf0ee')])
        s.configure('Accent.TButton',background=TEAL,foreground='white');s.map('Accent.TButton',background=[('active','#0c7d73'),('disabled','#aacdc5')])
        s.configure('TEntry',fieldbackground='white',padding=7,bordercolor=BORDER)
        s.configure('TCombobox',fieldbackground='white',padding=6)
        s.configure('Treeview',background=PAPER,fieldbackground=PAPER,rowheight=33,borderwidth=0,font=('DejaVu Sans',9));s.configure('Treeview.Heading',background='#e8efed',foreground=MUTED,font=('DejaVu Sans',9,'bold'),padding=8,relief='flat');s.map('Treeview',background=[('selected','#d3ece5')],foreground=[('selected',INK)])
        s.configure('Horizontal.TProgressbar',background=TEAL,troughcolor='#deebe5',borderwidth=0)
    def make_shell(self):
        sidebar=tk.Frame(self,bg=NAV,width=212);sidebar.pack(side='left',fill='y');sidebar.pack_propagate(False)
        logo=tk.Frame(sidebar,bg=NAV);logo.pack(fill='x',padx=22,pady=(30,8));tk.Label(logo,text='WES / LAB',font=('DejaVu Sans',20,'bold'),fg='white',bg=NAV).pack(anchor='w')
        tk.Label(sidebar,text='EXOME RESEARCH WORKSPACE',font=('DejaVu Sans',8),fg='#89aaa8',bg=NAV).pack(anchor='w',padx=23,pady=(0,28))
        self.nav={}
        for label in ('Projects','Files','Overview','Variants','Coverage','Reports','Run analysis'):
            b=tk.Button(sidebar,text=label,anchor='w',padx=23,pady=13,relief='flat',bd=0,highlightthickness=0,bg=NAV,fg='#adc5c4',activebackground='#24444b',activeforeground='white',font=('DejaVu Sans',11),cursor='hand2',command=lambda p=label:self.show_page(p));b.pack(fill='x',padx=10,pady=3);self.nav[label]=b
        tk.Frame(sidebar,bg='#355057',height=1).pack(fill='x',padx=22,pady=23)
        self.import_button=tk.Button(sidebar,text='+  Add analysis / VCF',bg=TEAL,fg='white',relief='flat',font=('DejaVu Sans',10,'bold'),padx=10,pady=12,cursor='hand2',command=self.import_dialog);self.import_button.pack(fill='x',padx=21,pady=(0,10))
        self.open_button=tk.Button(sidebar,text='Open saved project',bg='#24444b',fg='white',relief='flat',font=('DejaVu Sans',10),pady=11,cursor='hand2',command=self.choose_project);self.open_button.pack(fill='x',padx=21)
        self.example_button=tk.Button(sidebar,text='Open existing 571566',bg=NAV,fg='#9cc5bc',relief='flat',font=('DejaVu Sans',9),pady=12,command=self.open_example);self.example_button.pack(fill='x',padx=15,pady=7)
        tk.Label(sidebar,text='LOCAL BY DEFAULT\nFiles stay on this computer.',justify='left',bg=NAV,fg='#89aaa8',font=('DejaVu Sans',9)).pack(side='bottom',anchor='w',padx=23,pady=25)
        right=tk.Frame(self,bg=BG);right.pack(side='left',fill='both',expand=True)
        header=tk.Frame(right,bg=PAPER,highlightbackground=BORDER,highlightthickness=1);header.pack(fill='x')
        self.breadcrumb=tk.Label(header,text='Workspace / Overview',bg=PAPER,fg=MUTED,font=('DejaVu Sans',10),padx=26,pady=19);self.breadcrumb.pack(side='left')
        self.sample_label=tk.Label(header,text='No project open',bg=PAPER,fg=TEAL,font=('DejaVu Sans',10,'bold'));self.sample_label.pack(side='right',padx=24)
        bottom=tk.Frame(right,bg=PAPER,highlightbackground=BORDER,highlightthickness=1);bottom.pack(side='bottom',fill='x')
        self.status=tk.Label(bottom,text='Ready · research review, not a diagnostic result',anchor='w',bg=PAPER,fg=MUTED,font=('DejaVu Sans',9));self.status.pack(side='left',padx=20,pady=12)
        self.stop_button=ttk.Button(bottom,text='Stop task',command=self.cancel_task);self.stop_button.pack(side='right',padx=12,pady=5);self.stop_button.configure(state='disabled')
        self.progress=ttk.Progressbar(bottom,mode='indeterminate',length=120);self.progress.pack(side='right',padx=8)
        self.content=tk.Frame(right,bg=BG);self.content.pack(fill='both',expand=True,padx=25,pady=22)
    def title_block(self,parent,title,subtitle):
        tk.Label(parent,text=title,font=('DejaVu Sans',24,'bold'),bg=BG,fg=INK,anchor='w').pack(fill='x',pady=(0,7))
        tk.Label(parent,text=subtitle,font=('DejaVu Sans',10),bg=BG,fg=MUTED,anchor='w',justify='left',wraplength=950).pack(fill='x',pady=(0,20))
    def show_page(self,page):
        self.page=page;self.breadcrumb.configure(text='Workspace / '+page)
        for label,b in self.nav.items():b.configure(bg='#294a4e' if label==page else NAV,fg='white' if label==page else '#adc5c4')
        for widget in self.content.winfo_children():widget.destroy()
        if page=='Projects':self.people_page();return
        if page=='Files':self.files_page();return
        if page=='Run analysis':self.workflow_page();return
        if self.project is None:self.welcome();return
        if page=='Overview':self.overview()
        elif page=='Variants':self.variants_page()
        elif page=='Coverage':self.coverage_page()
        elif page=='Reports':self.reports_page()
    def welcome(self):
        self.title_block(self.content,'Your exome, made readable.','A private workspace for exploring WES variants, gene evidence and sequencing coverage.')
        card=tk.Frame(self.content,bg=PAPER,highlightbackground=BORDER,highlightthickness=1);card.pack(fill='x',pady=15)
        tk.Label(card,text='Start with a variant file',font=('DejaVu Sans',18,'bold'),bg=PAPER,fg=INK).pack(anchor='w',padx=30,pady=(30,12))
        tk.Label(card,text='Open .vcf or .vcf.gz files with SnpEff ANN or VEP CSQ annotations.\nUnannotated VCFs work for quality review; gene analysis needs annotations.',font=('DejaVu Sans',11),justify='left',bg=PAPER,fg=MUTED).pack(anchor='w',padx=30,pady=5)
        row=tk.Frame(card,bg=PAPER);row.pack(anchor='w',padx=30,pady=25)
        ttk.Button(row,text='Add WES / VCF',style='Accent.TButton',command=self.import_dialog).pack(side='left',padx=(0,10));ttk.Button(row,text='Open existing sample 571566',command=self.open_example).pack(side='left')
        tk.Label(self.content,text='Also available: an analysis-ready BAM workflow using the installed GRCh38 tools.\nFASTQ alignment is not included. SFARI membership is gene-level research evidence, not a diagnosis.',font=('DejaVu Sans',10),bg=BG,fg=MUTED,justify='left').pack(anchor='w',pady=22)
        if self.settings.get('recent'):
            tk.Label(self.content,text='RECENT PROJECTS',font=('DejaVu Sans',9,'bold'),bg=BG,fg=MUTED).pack(anchor='w',pady=(20,8))
            for p in self.settings['recent'][:4]:
                if Path(p).exists():ttk.Button(self.content,text=Path(p).name,command=lambda p=p:self.open_project(p)).pack(anchor='w',pady=4)
    def cards(self,parent,items):
        row=tk.Frame(parent,bg=BG);row.pack(fill='x',pady=(0,20))
        for i,(number,label,note) in enumerate(items):
            card=tk.Frame(row,bg=PAPER,highlightbackground=BORDER,highlightthickness=1);card.grid(row=0,column=i,sticky='nsew',padx=(0,12 if i<len(items)-1 else 0));row.columnconfigure(i,weight=1)
            tk.Label(card,text=label.upper(),bg=PAPER,fg=MUTED,font=('DejaVu Sans',8,'bold')).pack(anchor='w',padx=17,pady=(17,7));tk.Label(card,text=number,bg=PAPER,fg=INK,font=('DejaVu Sans',26,'bold')).pack(anchor='w',padx=17);tk.Label(card,text=note,bg=PAPER,fg=TEAL,font=('DejaVu Sans',8),wraplength=190,justify='left').pack(anchor='w',padx=17,pady=(6,17))
    def overview(self):
        page=ScrollPage(self.content);page.pack(fill='both',expand=True);p=page.body;s=self.summary_cache
        self.title_block(p,'A clearer view of your WES.',f"Sample {self.project.meta['sample']} · {self.project.meta['build']} · {self.project.meta['imported_at']}"+(f" · {self.project.meta['counts'].get('placeholder_alleles_skipped',0):,} deletion placeholders excluded" if self.project.meta['counts'].get('placeholder_alleles_skipped') else ''))
        self.cards(p,[(f"{s['total']:,}",'Variant alleles','Selected sample, imported records'),(f"{s['passing']:,}",'Passing calls','VCF FILTER = PASS'),(f"{s['sfari']:,}",'SFARI overlaps','Gene evidence, not pathogenicity'),(f"{s['clinical_assertions']:,}",'Clinical assertions','Require individual evidence review')])
        if self.project.meta['counts'].get('unannotated_alleles',0):
            n=self.project.meta['counts']['unannotated_alleles'];tk.Label(p,text=f'{n:,} alleles have no matched transcript annotation. Use Run analysis to annotate compatible GRCh38 VCFs, or supply ANN/CSQ data.',bg='#fff3df',fg=AMBER,font=('DejaVu Sans',10),wraplength=900,justify='left',padx=18,pady=14).pack(fill='x',pady=(0,18))
        charts=tk.Frame(p,bg=BG);charts.pack(fill='both',expand=True);charts.columnconfigure(0,weight=1);charts.columnconfigure(1,weight=1)
        impact={k:s['impact'].get(k,0) for k in ('HIGH','MODERATE','LOW','MODIFIER','UNKNOWN') if s['impact'].get(k,0)}
        Chart(charts,'Predicted transcript effects','Click a bar to inspect the corresponding variants.',impact,on_select=lambda x:self.explore(impact=x)).grid(row=0,column=0,sticky='nsew',padx=(0,12),pady=(0,14))
        depths={k:s['depths'].get(k,0) for k in ('<10','10–19','20–49','50–99','100+','Missing') if s['depths'].get(k,0)}
        Chart(charts,'Read depth at variant sites','Variant-site DP; this is not whole-exome coverage.',depths).grid(row=0,column=1,sticky='nsew',pady=(0,14))
        def chromkey(x):return (0,int(x)) if x.isdigit() else (1,x)
        chroms=dict(sorted(s['chromosomes'].items(),key=lambda x:chromkey(x[0])))
        Chart(charts,'Variants by chromosome','Scroll within this page to view all charts.',chroms,on_select=lambda x:self.explore(chromosome=x)).grid(row=1,column=0,sticky='nsew',padx=(0,12),pady=(0,14))
        Chart(charts,'Clinical database assertions','Conflicts remain separate from pathogenic assertions.',s['assertions']).grid(row=1,column=1,sticky='nsew',pady=(0,14))
        Chart(charts,'Variant types','Sequence differences in the imported VCF.',s['types']).grid(row=2,column=0,sticky='nsew',padx=(0,12),pady=(0,14))
        Chart(charts,'Measured gene coverage','Genes grouped by coding/flank bases reaching 20×.',s['coverage_bins']).grid(row=2,column=1,sticky='nsew',pady=(0,14))
        tk.Label(p,text='Representative effects prefer matching MANE transcripts when available. All transcript annotations remain available in variant details.\nDatabase assertions and predictions require allele, inheritance and phenotype review.',bg=BG,fg=MUTED,font=('DejaVu Sans',9),justify='left',wraplength=950).pack(anchor='w',pady=12)
    def explore(self,**filters):
        self.show_page('Variants')
        if filters.get('impact'):self.effect_var.set(filters['impact'])
        if filters.get('chromosome'):self.chrom_var.set(filters['chromosome'])
        if filters.get('query'):self.search_var.set(filters['query'])
        self.apply_filters()
    def variants_page(self):
        self.title_block(self.content,'Explore the variants.','Filter the evidence, inspect every transcript, and export exactly what you are reviewing.')
        tools=tk.Frame(self.content,bg=BG);tools.pack(fill='x',pady=(0,12))
        self.search_var=tk.StringVar();entry=ttk.Entry(tools,textvariable=self.search_var,width=29);entry.pack(side='left',padx=(0,8));entry.bind('<KeyRelease>',self.search_changed)
        tk.Label(tools,text='Gene / HGVS / locus',bg=BG,fg=MUTED,font=('DejaVu Sans',8)).pack(side='left',padx=(0,12))
        self.effect_var=tk.StringVar(value='All effects');effect=ttk.Combobox(tools,textvariable=self.effect_var,values=['All effects','HIGH','MODERATE','LOW','MODIFIER','UNKNOWN'],state='readonly',width=13);effect.pack(side='left',padx=4);effect.bind('<<ComboboxSelected>>',lambda _:self.apply_filters())
        self.chrom_var=tk.StringVar(value='All chromosomes');chrom=ttk.Combobox(tools,textvariable=self.chrom_var,values=['All chromosomes',*self.summary_cache['chromosomes']],state='readonly',width=16);chrom.pack(side='left',padx=4);chrom.bind('<<ComboboxSelected>>',lambda _:self.apply_filters())
        options=tk.Frame(self.content,bg=BG);options.pack(fill='x',pady=(0,12))
        self.pass_var=tk.BooleanVar(value=False);self.sfari_var=tk.BooleanVar(value=False);self.review_var=tk.BooleanVar(value=False)
        for text,var in [('PASS only',self.pass_var),('SFARI genes',self.sfari_var),('Review flags',self.review_var)]:ttk.Checkbutton(options,text=text,variable=var,command=self.apply_filters).pack(side='left',padx=(0,14))
        self.group_var=tk.StringVar(value='All review groups');groups=['All review groups','Clinical assertion','Predicted effect','Conflicting assertions','Transcript check','Support check','Filtered / unfiltered','Other'];group=ttk.Combobox(options,textvariable=self.group_var,values=groups,state='readonly',width=24);group.pack(side='left',padx=5);group.bind('<<ComboboxSelected>>',lambda _:self.apply_filters())
        ttk.Button(options,text='Export view',command=self.export_filtered).pack(side='right')
        ttk.Button(options,text='Inspect reads / IGV',command=self.open_igv_review).pack(side='right',padx=8)
        split=tk.PanedWindow(self.content,orient='vertical',bg=BORDER,sashwidth=5,relief='flat');split.pack(fill='both',expand=True)
        table_frame=tk.Frame(split,bg=PAPER);split.add(table_frame,minsize=180,stretch='always')
        self.table=self.tree(table_frame,[('gene','Gene(s)',155),('locus','Locus',140),('allele','REF → ALT',95),('effect','Effect',170),('gt','GT',55),('dp','DP',55),('gq','GQ',50),('clinical','ClinVar',160),('flag','Review group',145)])
        self.table.bind('<<TreeviewSelect>>',self.variant_selected);self.table.bind('<Double-1>',self.variant_window)
        detail_frame=tk.Frame(split,bg=PAPER);split.add(detail_frame,minsize=115,stretch='never')
        self.detail_text=ScrolledText(detail_frame,height=8,bg=PAPER,fg=INK,wrap='word',font=('DejaVu Sans',9),relief='flat',padx=15,pady=12);self.detail_text.pack(fill='both',expand=True);self.set_text(self.detail_text,'Select a variant to see genotype support, caveats and transcript context. Double-click for every transcript.')
        paging=tk.Frame(self.content,bg=BG);paging.pack(fill='x',pady=(10,0));self.row_count=tk.Label(paging,bg=BG,fg=MUTED,font=('DejaVu Sans',9));self.row_count.pack(side='left')
        self.next_button=ttk.Button(paging,text='Next →',command=lambda:self.change_page(1));self.next_button.pack(side='right');self.prev_button=ttk.Button(paging,text='← Previous',command=lambda:self.change_page(-1));self.prev_button.pack(side='right',padx=8)
        self.offset=0;self.apply_filters()
    def tree(self,parent,columns):
        frame=tk.Frame(parent,bg=PAPER);frame.pack(fill='both',expand=True)
        tree=ttk.Treeview(frame,columns=[c[0] for c in columns],show='headings',selectmode='browse');vertical=ttk.Scrollbar(frame,orient='vertical',command=tree.yview);horizontal=ttk.Scrollbar(frame,orient='horizontal',command=tree.xview);tree.configure(yscrollcommand=vertical.set,xscrollcommand=horizontal.set)
        tree.grid(row=0,column=0,sticky='nsew');vertical.grid(row=0,column=1,sticky='ns');horizontal.grid(row=1,column=0,sticky='ew');frame.columnconfigure(0,weight=1);frame.rowconfigure(0,weight=1)
        for key,label,width in columns:tree.heading(key,text=label);tree.column(key,width=width,minwidth=40,stretch=True)
        tree.tag_configure('assertion',foreground='#a76021');tree.tag_configure('muted',foreground=MUTED);return tree
    def search_changed(self,_=None):
        if self.debounce:self.after_cancel(self.debounce)
        self.debounce=self.after(250,self.apply_filters)
    def filters(self):
        return dict(query=self.search_var.get(),pass_only=self.pass_var.get(),sfari_only=self.sfari_var.get(),review_only=self.review_var.get(),impact=self.effect_var.get(),chromosome=self.chrom_var.get(),group=self.group_var.get())
    def apply_filters(self,reset=True):
        self.debounce=None
        if self.page!='Variants' or self.project is None:return
        if reset:self.offset=0
        count,rows=self.project.variants(self.filters(),offset=self.offset,limit=self.page_size,sort='review' if self.review_var.get() else 'position');self.visible_rows={str(r['id']):r for r in rows};self.table.delete(*self.table.get_children())
        for r in rows:
            allele=f"{r['ref']} → {r['alt']}";allele=allele if len(allele)<32 else allele[:29]+'…'
            values=(r['genes'] or '—',f"{r['chromosome']}:{r['position']:,}",allele,r['effect'],r['genotype'],self.num(r['depth']),self.num(r['gq']),r['clinvar'].replace('_',' ') or 'Not annotated',r['review_group'])
            self.table.insert('','end',iid=str(r['id']),values=values,tags=('assertion',) if r['review_group']=='Clinical assertion' else ())
        self.row_count.configure(text=f"{count:,} matches · showing {self.offset+1 if rows else 0:,}–{self.offset+len(rows):,}");self.prev_button.configure(state='normal' if self.offset else 'disabled');self.next_button.configure(state='normal' if self.offset+len(rows)<count else 'disabled')
    def change_page(self,direction):self.offset=max(0,self.offset+direction*self.page_size);self.apply_filters(reset=False)
    @staticmethod
    def num(value):return '—' if value is None else str(int(value)) if float(value).is_integer() else f'{value:.2f}'
    @staticmethod
    def set_text(widget,text):widget.configure(state='normal');widget.delete('1.0','end');widget.insert('1.0',text);widget.configure(state='disabled')
    def variant_description(self,r,full=False):
        af='Unavailable (not proof of rarity)' if r['population_af'] is None else f"{100*r['population_af']:.5g}%"
        text=f"{r['genes'] or 'No gene annotation'}  /  {r['chromosome']}:{r['position']} {r['ref']}>{r['alt']}\n{r['hgvs_c']}  {r['hgvs_p']}  ·  {r['transcript']}\nGT {r['genotype']}   DP {self.num(r['depth'])}   GQ {self.num(r['gq'])}   AD {r['allele_depths']}   FILTER {r['filter']}\nClinVar: {r['clinvar'] or 'Not annotated'}\nReview status: {r['clinvar_review'] or 'Not supplied'}\nConditions attached to assertions: {r['condition'].replace('_',' ') or 'Not supplied'}\nPopulation ALT frequency: {af}\nReview notes: {r['notes'] or 'No automated notes; this is not confirmation.'}"
        if full:
            text+='\n\nALL TRANSCRIPT ANNOTATIONS\n'
            for a in r['annotations']:text+='\n'+('MANE · ' if a.get('mane') else '')+f"{a['gene']}  {a['transcript']}\n  {a['effect']} ({a['impact']}) · {a['hgvs_c']} · {a['hgvs_p']}\n"+(f"  Notes: {a['warnings']}\n" if a.get('warnings') else '')
        return text
    def variant_selected(self,_=None):
        ids=self.table.selection()
        if ids:
            r=self.project.detail(int(ids[0]));self.set_text(self.detail_text,self.variant_description(r))
    def open_igv_review(self):
        import subprocess,sys
        source=Path(self.project.meta.get('source','')).resolve()
        if source!=(ROOT/'results/wgs/571566.annotated.vcf.gz').resolve():
            messagebox.showinfo('IGV bundle','This configured read-review bundle belongs to sample 571566. Open that analysis to inspect its reads.');return
        ids=self.table.selection();r=self.project.detail(int(ids[0])) if ids else None
        args=[sys.executable,'-m','wes_lab.igv_review']
        if r:args+=['--chrom',str(r['chromosome']),'--position',str(r['position'])]
        subprocess.Popen(args,cwd=ROOT)
    def variant_window(self,_=None):
        ids=self.table.selection()
        if not ids:return
        r=self.project.detail(int(ids[0]));win=tk.Toplevel(self);win.title(f"Variant details / {r['genes'] or r['chromosome']}");win.geometry('900x650');win.configure(bg=BG)
        text=ScrolledText(win,wrap='word',font=('DejaVu Sans',10),padx=20,pady=20,bg=PAPER,fg=INK);text.pack(fill='both',expand=True,padx=12,pady=12);self.set_text(text,self.variant_description(r,True))
        row=tk.Frame(win,bg=BG);row.pack(fill='x',padx=15,pady=(0,15));ttk.Button(row,text='Copy allele',command=lambda:(self.clipboard_clear(),self.clipboard_append(f"{r['chromosome']}:{r['position']} {r['ref']}>{r['alt']}"))).pack(side='left')
        term=(r['genes'].split(',')[0]+' '+r['hgvs_c']).strip()
        ttk.Button(row,text='Search public ClinVar',command=lambda:webbrowser.open('https://www.ncbi.nlm.nih.gov/clinvar/?term='+quote(term))).pack(side='left',padx=10)
    def coverage_page(self):
        self.title_block(self.content,'Know what was covered.','Coverage comes from a measured table—not from the presence or absence of variant calls.')
        summary=self.summary_cache;n=summary['coverage_genes'];self.cards(self.content,[(f'{n:,}','Genes measured','From the attached coverage table'),(f"{summary['coverage_complete']:,}",'At least 90% at 20×','Coding/flank interval completeness'),(f"{n-summary['coverage_complete']:,}",'Coverage gaps','Below 90% of supplied bases at 20×')])
        controls=tk.Frame(self.content,bg=BG);controls.pack(fill='x',pady=(0,14));self.coverage_query=tk.StringVar();ttk.Entry(controls,textvariable=self.coverage_query,width=25).pack(side='left');self.coverage_query.trace_add('write',lambda *_:self.refresh_coverage());self.gaps_var=tk.BooleanVar(value=False);ttk.Checkbutton(controls,text='Coverage gaps only',variable=self.gaps_var,command=self.refresh_coverage).pack(side='left',padx=14);ttk.Button(controls,text='Attach coverage TSV',command=self.attach_coverage).pack(side='right');ttk.Button(controls,text='Export table',command=self.export_coverage).pack(side='right',padx=9)
        self.coverage_tree=self.tree(self.content,[('gene','Gene',160),('bases','Measured bases',150),('mean','Mean depth',140),('pct','Bases at 20× (%)',170),('source','Table source',260)]);self.coverage_tree.bind('<Double-1>',self.coverage_gene)
        self.coverage_label=tk.Label(self.content,bg=BG,fg=MUTED,font=('DejaVu Sans',9),anchor='w');self.coverage_label.pack(fill='x',pady=10)
        self.refresh_coverage()
    def refresh_coverage(self):
        if self.page!='Coverage' or not hasattr(self,'coverage_tree'):return
        rows=self.project.coverage(self.coverage_query.get(),self.gaps_var.get());self.coverage_tree.delete(*self.coverage_tree.get_children())
        for r in rows[:600]:self.coverage_tree.insert('','end',iid=r['gene'],values=(r['gene'],r['bases'] or '—',f"{r['mean_depth']:.2f}×",f"{r['percent_ge20']:.2f}%",r['source']),tags=('assertion',) if r['percent_ge20']<90 else ())
        self.coverage_label.configure(text=f'{len(rows):,} matches · showing first {min(600,len(rows)):,}. Double-click a gene to explore its variants.' if rows else 'No measured coverage attached. Import a TSV with gene, mean_depth and percent_ge20 columns.')
    def coverage_gene(self,_=None):
        ids=self.coverage_tree.selection()
        if ids:self.explore(query=ids[0])
    def attach_coverage(self):
        if self.busy:return
        path=filedialog.askopenfilename(title='Attach measured gene coverage',filetypes=[('TSV tables','*.tsv'),('All files','*')])
        if not path:return
        try:self.project.import_coverage(path,source=Path(path).name);self.summary_cache=self.project.summary();self.show_page('Coverage');self.status.configure(text='Measured coverage attached to this project.')
        except Exception as e:messagebox.showerror('Coverage import failed',str(e))
    def export_coverage(self):
        path=filedialog.asksaveasfilename(title='Export measured coverage',initialdir=self.person.path/'reports' if self.person else self.workspace,defaultextension='.tsv',filetypes=[('TSV table','*.tsv')])
        if not path:return
        try:
            rows=self.project.coverage()
            with Path(path).open('w',newline='',encoding='utf-8') as f:
                w=csv.DictWriter(f,['gene','bases','mean_depth','percent_ge20','source'],delimiter='\t');w.writeheader();w.writerows(rows)
            self.person.register_file(path,kind='Coverage table') if self.person else None;self.status.configure(text='Coverage table exported.')
        except Exception as e:messagebox.showerror('Export failed',str(e))
    def reports_page(self):
        self.title_block(self.content,'Take the evidence with you.','Export a readable report, share a filtered table, or revisit the source analysis.')
        card=tk.Frame(self.content,bg=PAPER,highlightbackground=BORDER,highlightthickness=1);card.pack(fill='x',pady=5)
        tk.Label(card,text='Project report',bg=PAPER,fg=INK,font=('DejaVu Sans',16,'bold')).pack(anchor='w',padx=25,pady=(25,10));tk.Label(card,text='A standalone HTML report with graphs, counts, coverage summary and database-assertion candidates.\nOpens in a browser and can be printed to PDF. It contains genomic information; share deliberately.',bg=PAPER,fg=MUTED,font=('DejaVu Sans',10),justify='left',wraplength=850).pack(anchor='w',padx=25,pady=6)
        ttk.Button(card,text='Export HTML report',style='Accent.TButton',command=self.export_report).pack(anchor='w',padx=25,pady=22)
        row=tk.Frame(self.content,bg=BG);row.pack(fill='x',pady=20);ttk.Button(row,text='Export all variants as CSV',command=lambda:self.export_filtered(all_rows=True)).pack(side='left');ttk.Button(row,text='Open project metadata',command=lambda:self.text_file(self.project.path/'project.json')).pack(side='left',padx=12)
        if self.is_existing_sample():
            tk.Label(self.content,text='EXISTING RESEARCH REVIEWS / 571566',bg=BG,fg=MUTED,font=('DejaVu Sans',9,'bold')).pack(anchor='w',pady=(15,10))
            for title,path in [('Expanded exome review — FDXR, WNT10A and other findings',ROOT/'results/wgs/exome_review/REVIEW.md'),('Clinical handoff',ROOT/'results/wgs/exome_review/CLINICAL_HANDOFF.md'),('Original 18-candidate review',ROOT/'results/wgs/candidate_review/REVIEW.md')]:
                if path.exists():ttk.Button(self.content,text=title,command=lambda p=path:self.text_file(p)).pack(anchor='w',pady=5)
        tk.Label(self.content,text='Interpretation limits\n\nA VCF import does not call copy-number changes, repeat expansions, structural variants or mosaicism.\nSFARI is gene-level research evidence. A rare or HIGH-effect variant is not automatically pathogenic.\nNo “autism probability” or diagnostic score is calculated.',bg=BG,fg=MUTED,font=('DejaVu Sans',10),justify='left',wraplength=900).pack(anchor='w',pady=25)
    def export_filtered(self,all_rows=False):
        if not self.project or self.busy:return
        path=filedialog.asksaveasfilename(title='Export variants',initialdir=self.person.path/'reports' if self.person else self.workspace,defaultextension='.csv',filetypes=[('CSV table','*.csv')])
        if not path:return
        filters={} if all_rows or self.page!='Variants' else self.filters()
        owner=self.person
        def done(_):
            if owner:owner.register_file(path,kind='Variant export')
            self.status.configure(text='Variant table exported.')
        self.start_task(lambda:self.project.export_csv(path,filters),done,'Exporting variants…')
    def export_report(self):
        if not self.project or self.busy:return
        path=filedialog.asksaveasfilename(title='Export project report',initialdir=self.person.path/'reports' if self.person else self.workspace,defaultextension='.html',filetypes=[('HTML report','*.html')])
        if path:
            try:
                self.project.export_report(path)
                if self.person:self.person.register_file(path,kind='Report')
                self.status.configure(text=f'Report saved: {path}')
            except Exception as e:messagebox.showerror('Report export failed',str(e))
    def text_file(self,path):
        win=tk.Toplevel(self);win.title(Path(path).name);win.geometry('950x740');text=ScrolledText(win,wrap='word',font=('DejaVu Sans',10),bg=PAPER,fg=INK,padx=20,pady=20);text.pack(fill='both',expand=True)
        try:self.set_text(text,Path(path).read_text())
        except OSError as e:self.set_text(text,str(e))
    def is_existing_sample(self):return self.project and Path(self.project.meta['source'])==ROOT/'results/wgs/571566.annotated.vcf.gz'
    def import_dialog(self,path=None):
        if self.busy:return
        if not self.person:self.new_person_dialog(on_created=lambda:self.import_dialog(path));return
        path=path or filedialog.askopenfilename(title='Open WES variant file',filetypes=[('VCF files','*.vcf *.vcf.gz'),('All files','*')])
        if not path:return
        try:header=read_header(path)
        except Exception as e:messagebox.showerror('Cannot open VCF',str(e));return
        dialog=tk.Toplevel(self);dialog.title('Import WES analysis');dialog.geometry('570x390');dialog.configure(bg=BG);dialog.transient(self);dialog.grab_set()
        tk.Label(dialog,text='Add a WES analysis',font=('DejaVu Sans',17,'bold'),bg=BG,fg=INK).pack(anchor='w',padx=25,pady=(25,15));tk.Label(dialog,text=Path(path).name,bg=BG,fg=MUTED,wraplength=500).pack(anchor='w',padx=25)
        form=tk.Frame(dialog,bg=BG);form.pack(fill='x',padx=25,pady=15)
        sample=tk.StringVar(value=header['samples'][0] if header['samples'] else 'Sites only');build=tk.StringVar(value='Unknown');name=tk.StringVar(value=Path(path).name.split('.vcf')[0])
        for i,(label,var,values) in enumerate([('Analysis name',name,None),('Sample',sample,header['samples'] or ['Sites only']),('Genome build',build,['Unknown','GRCh38','GRCh37'])]):
            tk.Label(form,text=label,bg=BG,fg=MUTED).grid(row=i,column=0,sticky='w',pady=9)
            widget=ttk.Combobox(form,textvariable=var,values=values,state='readonly') if values else ttk.Entry(form,textvariable=var)
            widget.grid(row=i,column=1,sticky='ew',padx=(15,0));form.columnconfigure(1,weight=1)
        tk.Label(dialog,text='Use the actual reference build. GRCh38 MANE comparisons are applied only when GRCh38 is selected.',bg=BG,fg=MUTED,wraplength=500,font=('DejaVu Sans',9),justify='left').pack(anchor='w',padx=25)
        def confirm():
            safe=re_safe(name.get())
            if not safe:messagebox.showerror('Analysis name required','Enter an analysis name.',parent=dialog);return
            target=self.person.new_analysis_directory(safe)
            selected_sample=sample.get() if header['samples'] else None;selected_build=build.get();dialog.destroy();self.begin_import(path,target,selected_sample,selected_build)
        ttk.Button(dialog,text='Import variants',style='Accent.TButton',command=confirm).pack(anchor='e',padx=25,pady=18)
    def begin_import(self,source,target,sample=None,build='Unknown',example=False):
        owner=self.person
        def job():
            result=import_vcf(source,target,sample,build,progress=lambda x:self.events.put(('progress',x)),cancel=self.cancel)
            if example:
                cov=ROOT/'results/wgs/exome_review/all_gene_coverage.tsv'
                if not cov.exists():cov=ROOT/'results/wgs/qc/571566.sfari_gene_coverage.tsv'
                if cov.exists():result.import_coverage(cov,source=cov.name)
                # Supplement missing MANE genes from the existing SFARI/Ensembl table.
                supplement=ROOT/'results/wgs/qc/571566.sfari_gene_coverage.tsv'
                if supplement.exists() and cov!=supplement:
                    with result.connect() as db,supplement.open() as f:
                        for r in csv.DictReader(f,delimiter='\t'):
                            db.execute('INSERT OR IGNORE INTO coverage VALUES (?,?,?,?,?)',(r['gene'],int(r['coding_plus_flank_bases']),float(r['mean_depth']),float(r['percent_ge20']),'SFARI MANE / Ensembl supplement'))
            return result
        self.start_task(job,lambda project:self.loaded(project,person=owner),'Importing variants into a local project…')
    def open_example(self):
        if self.busy:return
        source=ROOT/'results/wgs/571566.annotated.vcf.gz';target=self.workspace/'571566'
        if not source.exists():messagebox.showinfo('Existing sample unavailable','The bundled sample path is unavailable. Add your own VCF.');return
        if (target/'variants.sqlite').exists():self.open_project(target)
        else:
            owner=find_analysis_owner(self.people_root,target) or PersonProject.create(self.people_root,'571566')
            self.open_person(owner);self.begin_import(source,target,'571566','GRCh38',example=True)
    def choose_project(self):
        if self.busy:return
        path=filedialog.askdirectory(title='Open person project or saved analysis',initialdir=self.people_root)
        if path:self.open_project(path)
    def open_project(self,path):
        if self.busy:return
        if (Path(path)/'person.json').exists():self.open_person(path);return
        try:self.loaded(Project(path))
        except Exception as e:messagebox.showerror('Could not open project',str(e))
    def loaded(self,project,person=None):
        person=person or find_analysis_owner(self.people_root,project.path)
        if person is None:person=PersonProject.create(self.people_root,str(project.meta['sample'])[:150] or 'Imported person')
        person.register_analysis(project.path,project.meta['sample'],project.meta['build'],label=project.path.name)
        source=Path(project.meta['source'])
        if source.exists():person.register_file(source,source='Source VCF for analysis')
        self.person=person;self.project=project;self.summary_cache=project.summary();self.sample_label.configure(text=f"{person.data['name']} / {project.meta['sample']}")
        if self.is_existing_sample():
            for path in (ROOT/'results/wgs/571566.recal.bam',ROOT/'results/wgs/571566.recal.bam.bai',ROOT/'results/wgs/571566.recal.bai',ROOT/'results/wgs/571566.annotated.vcf.gz.tbi',ROOT/'results/wgs/exome_review/REVIEW.md',ROOT/'results/wgs/exome_review/CLINICAL_HANDOFF.md'):
                if path.exists():person.link_file(path)
        self.remember_person(person)
        path=str(project.path);self.settings['recent']=[path]+[p for p in self.settings.get('recent',[]) if p!=path];self.settings_path.write_text(json.dumps(self.settings,indent=2)+'\n');self.show_page('Overview');self.status.configure(text=f"{person.data['name']} · {self.summary_cache['total']:,} variant alleles indexed locally")
    def workflow_page(self):
        if not self.person:
            self.title_block(self.content,'Choose a person project.','Analysis inputs and outputs belong to that person.');ttk.Button(self.content,text='Create person project',style='Accent.TButton',command=self.new_person_dialog).pack(anchor='w');return
        self.title_block(self.content,'Run a local analysis.','Sequential jobs with saved sections, live progress and a restartable workflow folder.')
        card=tk.Frame(self.content,bg=PAPER,highlightbackground=BORDER,highlightthickness=1);card.pack(fill='x',pady=5)
        tk.Label(card,text='Installed GRCh38 workflow',bg=PAPER,fg=INK,font=('DejaVu Sans',16,'bold')).pack(anchor='w',padx=24,pady=(23,12))
        missing=probe_bundle();available='Ready — installed local resources found' if not missing else 'Unavailable — missing: '+', '.join(missing)
        tk.Label(card,text=available,bg=PAPER,fg=TEAL if not missing else AMBER,font=('DejaVu Sans',10),wraplength=850).pack(anchor='w',padx=24,pady=5)
        tk.Label(card,text='BAM: supply an analysis-ready, coordinate-sorted, indexed BAM with one sample.\nAlignment, duplicate marking and recalibration assessment happen before this workflow.\nVCF: annotate a single-sample, accession-named GRCh38 VCF.\nTargets: installed Ensembl exons +100bp, an approximation rather than a capture-kit design.\nResources: one job at a time; calling heap 2 GB, annotation heap 3 GB. X/Y calls remain provisional.',bg=PAPER,fg=MUTED,font=('DejaVu Sans',10),justify='left',wraplength=880).pack(anchor='w',padx=24,pady=14)
        form=tk.Frame(card,bg=PAPER);form.pack(fill='x',padx=24,pady=10);self.workflow_mode=tk.StringVar(value='Analysis-ready BAM');self.workflow_input=tk.StringVar();self.workflow_output=tk.StringVar(value=str(self.person.path/'analyses'/('workflow_'+time.strftime('%Y%m%d_%H%M%S'))))
        tk.Label(form,text='Input type',bg=PAPER,fg=MUTED).grid(row=0,column=0,sticky='w',pady=10);ttk.Combobox(form,textvariable=self.workflow_mode,values=['Analysis-ready BAM','VCF annotation'],state='readonly',width=24).grid(row=0,column=1,sticky='w',padx=12)
        for row,(label,var) in enumerate([('Input file',self.workflow_input),('Workflow folder',self.workflow_output)],1):
            tk.Label(form,text=label,bg=PAPER,fg=MUTED).grid(row=row,column=0,sticky='w',pady=10);ttk.Entry(form,textvariable=var,width=55).grid(row=row,column=1,sticky='ew',padx=12);ttk.Button(form,text='Browse',command=lambda var=var,row=row:self.workflow_browse(var,row)).grid(row=row,column=2)
        form.columnconfigure(1,weight=1)
        controls=tk.Frame(card,bg=PAPER);controls.pack(fill='x',padx=24,pady=(10,23));self.workflow_start=ttk.Button(controls,text='Start / resume workflow',style='Accent.TButton',command=self.run_workflow);self.workflow_start.pack(side='left');ttk.Button(controls,text='View workflow log',command=self.view_workflow_log).pack(side='left',padx=10)
        if missing or self.busy:self.workflow_start.configure(state='disabled')
        tk.Label(self.content,text='Resume by selecting the same input and workflow folder. Completed sections are reused only if the configuration matches.\nImported files and completed outputs are retained when you stop a task. FASTQ alignment and clinical interpretation are not automated here.',bg=BG,fg=MUTED,font=('DejaVu Sans',9),wraplength=950,justify='left').pack(anchor='w',pady=20)
    def workflow_browse(self,var,row):
        if self.busy:return
        path=filedialog.askdirectory(title='Choose workflow folder') if row==2 else filedialog.askopenfilename(title='Select analysis input',filetypes=[('Analysis inputs','*.bam *.vcf *.vcf.gz'),('All files','*')])
        if path:var.set(path)
    def run_workflow(self):
        if self.busy:return
        source=self.workflow_input.get().strip();target=self.workflow_output.get().strip();mode='bam' if self.workflow_mode.get()=='Analysis-ready BAM' else 'vcf'
        if not source or not target:messagebox.showerror('Input needed','Choose an input file and workflow folder.');return
        owner=self.person
        if not Path(target).resolve().is_relative_to(owner.path):messagebox.showerror('Choose this person’s folder','Keep workflow outputs under this person project to avoid mixing samples.');return
        self.workflow_start.configure(state='disabled')
        def run():return Pipeline(source,target,mode,progress=lambda x:self.events.put(('progress',x)),cancel=self.cancel).run()
        def done(path):
            output=Path(target)/'viewer_project';self.begin_import(path,output,build='GRCh38')
        self.start_task(run,done,'Preparing the local analysis workflow…')
    def view_workflow_log(self):
        p=Path(self.workflow_output.get())/'workflow.log'
        if p.exists():self.text_file(p)
        else:messagebox.showinfo('No log yet','A log appears after the workflow starts.')
    def start_task(self,function,on_success,message):
        if self.busy:return
        self.busy=True;self.cancel=threading.Event();self.on_success=on_success;self.status.configure(text=message);self.progress.start(12);self.stop_button.configure(state='normal');self.import_button.configure(state='disabled');self.open_button.configure(state='disabled');self.example_button.configure(state='disabled')
        def worker():
            try:self.events.put(('done',function()))
            except Cancelled as e:self.events.put(('cancelled',str(e) or 'Stopped; completed sections retained.'))
            except Exception as e:self.events.put(('error',(str(e),traceback.format_exc())))
        self.worker=threading.Thread(target=worker,daemon=True);self.worker.start()
    def poll(self):
        try:
            while True:
                event,value=self.events.get_nowait()
                if event=='progress':self.status.configure(text=str(value)[:125])
                else:
                    self.busy=False;self.progress.stop();self.stop_button.configure(state='disabled');self.import_button.configure(state='normal');self.open_button.configure(state='normal');self.example_button.configure(state='normal')
                    if event=='done':self.on_success(value)
                    elif event=='error':
                        self.status.configure(text='Task failed · source files retained');(self.workspace/'last_error.log').write_text(value[1]);messagebox.showerror('Task could not finish',value[0])
                    else:self.status.configure(text=value)
                    if self.page=='Run analysis' and hasattr(self,'workflow_start') and self.workflow_start.winfo_exists():self.workflow_start.configure(state='disabled' if self.busy or probe_bundle() else 'normal')
        except queue.Empty:pass
        self.after(100,self.poll)
    def cancel_task(self):
        if self.busy:self.cancel.set();self.status.configure(text='Stopping after a safe checkpoint…');self.stop_button.configure(state='disabled')
    def close(self):
        if self.busy:
            if messagebox.askyesno('Task in progress','Stop the current task and close? Completed workflow sections will be saved.'):
                self.cancel_task();self.wait_to_close()
        else:self.destroy()
    def wait_to_close(self):
        if self.worker and self.worker.is_alive():self.after(200,self.wait_to_close)
        else:self.destroy()

def re_safe(text):
    import re
    return re.sub(r'[^A-Za-z0-9._-]+','_',text.strip()).strip('._')[:80]

def main():
    parser=argparse.ArgumentParser(description='WES Lab desktop research viewer');parser.add_argument('--project',type=Path);parser.add_argument('--example',action='store_true');args=parser.parse_args()
    app=App(example=args.example,project=args.project);app.mainloop()

if __name__=='__main__':main()
