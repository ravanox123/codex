"""Person-project and file-management screens for the native desktop app."""
import json,threading,subprocess
from pathlib import Path
import tkinter as tk
from tkinter import ttk,filedialog,messagebox
from tkinter.scrolledtext import ScrolledText
from .people import PersonProject,list_people,find_analysis_owner,file_kind,public_url
from .data import Project,Cancelled
BG='#f3f6f5';PAPER='#ffffff';INK='#18373e';MUTED='#637c80';TEAL='#158f83';BORDER='#dce6e3'

def size_text(value):
    value=float(value)
    for unit in ('B','KB','MB','GB','TB'):
        if value<1024 or unit=='TB':return f'{value:,.1f} {unit}'
        value/=1024

class PeopleUI:
    def people_page(self):
        self.title_block(self.content,'One person. One project.','Keep each person’s FASTQ, BAM, variants, analyses and reports together.')
        row=tk.Frame(self.content,bg=BG);row.pack(fill='x',pady=(0,18))
        ttk.Button(row,text='+ New person project',style='Accent.TButton',command=self.new_person_dialog).pack(side='left')
        ttk.Button(row,text='Open project folder',command=self.choose_project).pack(side='left',padx=10)
        self.people_table=self.tree(self.content,[('name','Person / identifier',230),('files','Files',75),('analyses','Analyses',90),('updated','Last updated',180),('path','Project folder',300)])
        self.people_entries={}
        for person in list_people(self.people_root):
            d=person.data;self.people_entries[d['id']]=person
            self.people_table.insert('','end',iid=d['id'],values=(d['name'],len(d['files']),len(d['analyses']),d['updated_at'],str(person.path)))
        self.people_table.bind('<Double-1>',lambda _:self.open_selected_person())
        ttk.Button(self.content,text='Open selected person',command=self.open_selected_person).pack(anchor='e',pady=14)
        tk.Label(self.content,text='Use a name or anonymous identifier. New downloads and analyses are stored inside the selected person’s folder.\nLarge files already on your computer can be linked in place or copied into the project.',bg=BG,fg=MUTED,font=('DejaVu Sans',10),justify='left',wraplength=950).pack(anchor='w',pady=12)
    def open_selected_person(self):
        ids=self.people_table.selection()
        if ids:self.open_person(self.people_entries[ids[0]])
    def new_person_dialog(self,on_created=None):
        if self.busy:return
        self.edit_person_dialog(None,on_created)
    def edit_person_dialog(self,person=None,on_created=None):
        if self.busy:return
        win=tk.Toplevel(self);win.title('Person project');win.geometry('580x380');win.configure(bg=BG);win.transient(self);win.grab_set()
        tk.Label(win,text='Person project' if person else 'Create a person project',bg=BG,fg=INK,font=('DejaVu Sans',18,'bold')).pack(anchor='w',padx=25,pady=(23,12))
        tk.Label(win,text='Person name or anonymous identifier',bg=BG,fg=MUTED).pack(anchor='w',padx=25)
        name=tk.StringVar(value=person.data['name'] if person else '');ttk.Entry(win,textvariable=name).pack(fill='x',padx=25,pady=9)
        tk.Label(win,text='Notes (optional)',bg=BG,fg=MUTED).pack(anchor='w',padx=25,pady=(8,0))
        notes=ScrolledText(win,height=5,font=('DejaVu Sans',10),wrap='word');notes.pack(fill='x',padx=25,pady=8)
        if person:notes.insert('1.0',person.data.get('notes',''))
        def save():
            try:
                if person:person.update(name.get(),notes.get('1.0','end').strip());result=person
                else:result=PersonProject.create(self.people_root,name.get(),notes.get('1.0','end').strip())
                win.destroy();self.open_person(result)
                if on_created:on_created()
            except Exception as e:messagebox.showerror('Could not save person',str(e),parent=win)
        ttk.Button(win,text='Save project' if person else 'Create project',style='Accent.TButton',command=save).pack(anchor='e',padx=25,pady=16)
    def remember_person(self,person):
        self.settings['last_person']=str(person.path)
        self.settings_path.write_text(json.dumps(self.settings,indent=2)+'\n')
    def open_person(self,person,page='Files'):
        if self.busy:return
        try:
            person=person if isinstance(person,PersonProject) else PersonProject(person)
            self.person=person;self.project=None;self.summary_cache={};self.sample_label.configure(text=person.data['name']);self.remember_person(person)
            analyses=person.read()['analyses'];active=person.data.get('active_analysis');selected=next((a for a in analyses if a['id']==active),analyses[-1] if analyses else None)
            if selected and (Path(selected['path'])/'variants.sqlite').exists():self.loaded(Project(selected['path']),person=person)
            self.show_page(page);self.status.configure(text='Person project: '+person.data['name'])
        except Exception as e:messagebox.showerror('Could not open person project',str(e))
    def files_page(self):
        if not self.person:
            self.title_block(self.content,'Choose a person first.','Files and downloads belong to one person project.')
            ttk.Button(self.content,text='Create person project',style='Accent.TButton',command=self.new_person_dialog).pack(anchor='w')
            ttk.Button(self.content,text='Browse projects',command=lambda:self.show_page('Projects')).pack(anchor='w',pady=12);return
        self.person.data=self.person.read();d=self.person.data
        self.title_block(self.content,d['name']+' / files','FASTQ, alignments, variants and reports, organized under this person’s project.')
        actions=tk.Frame(self.content,bg=BG);actions.pack(fill='x',pady=(0,14))
        for label,command,accent in [('Download links',self.download_dialog,True),('Link existing files',lambda:self.add_person_files(False),False),('Copy files here',lambda:self.add_person_files(True),False),('Open folder',self.open_person_folder,False),('Edit person',lambda:self.edit_person_dialog(self.person),False)]:
            ttk.Button(actions,text=label,command=command,style='Accent.TButton' if accent else 'TButton').pack(side='left',padx=(0,8))
        self.files_tabs=ttk.Notebook(self.content);self.files_tabs.pack(fill='both',expand=True)
        files_tab=tk.Frame(self.files_tabs,bg=BG);downloads_tab=tk.Frame(self.files_tabs,bg=BG);analyses_tab=tk.Frame(self.files_tabs,bg=BG)
        self.files_tabs.add(files_tab,text=f"Files ({len(d['files'])})");self.files_tabs.add(downloads_tab,text='Downloads');self.files_tabs.add(analyses_tab,text=f"Analyses ({len(d['analyses'])})")
        self.file_table=self.tree(files_tab,[('name','Filename',290),('kind','Type',120),('storage','Storage',100),('size','Size',100),('status','State',140)])
        self.file_entries={}
        for r in self.person.files():
            self.file_entries[r['id']]=r;state='Missing' if not r['exists'] else 'Changed since added' if r['changed'] else 'Ready'
            self.file_table.insert('','end',iid=r['id'],values=(r['name'],r['kind'],r['storage'],size_text(r['size']),state))
        self.file_note=tk.Label(files_tab,text='Select a file for its location and checksum. Linked files remain in their original location.',bg=BG,fg=MUTED,justify='left',anchor='w',wraplength=930,font=('DejaVu Sans',9));self.file_note.pack(fill='x',pady=12)
        self.file_table.bind('<<TreeviewSelect>>',self.person_file_selected)
        row=tk.Frame(files_tab,bg=BG);row.pack(fill='x',pady=(0,12));ttk.Button(row,text='Use selected file in analysis',command=self.use_person_file).pack(side='left');ttk.Button(row,text='Copy file path',command=self.copy_person_path).pack(side='left',padx=10);ttk.Button(row,text='Refresh files',command=lambda:self.show_page('Files')).pack(side='right')
        self.download_table=self.tree(downloads_tab,[('name','Filename',250),('status','Status',120),('progress','Downloaded',160),('source','Source',350)])
        self.download_entries={}
        for job in self.person.downloads():
            self.download_entries[job['id']]=job;progress=size_text(job['bytes'])+(' / '+size_text(job['total']) if job.get('total') else '')
            self.download_table.insert('','end',iid=job['id'],values=(job['filename'],job['status'],progress,job['source_display']))
        self.download_note=tk.Label(downloads_tab,text='Select a download to inspect its status or error.',bg=BG,fg=MUTED,font=('DejaVu Sans',9),anchor='w',justify='left',wraplength=930);self.download_note.pack(fill='x',pady=8)
        def download_selected(_=None):
            ids=self.download_table.selection()
            if ids:
                job=self.download_entries[ids[0]];self.download_note.configure(text=(job.get('error') or ('State: '+job['status'])).replace(job['url'],job['source_display']))
        self.download_table.bind('<<TreeviewSelect>>',download_selected)
        row=tk.Frame(downloads_tab,bg=BG);row.pack(fill='x',pady=12);ttk.Button(row,text='Resume selected download',command=self.resume_selected_download).pack(side='left');ttk.Button(row,text='Resume all pending',command=self.resume_pending_downloads).pack(side='left',padx=10);ttk.Button(row,text='Refresh',command=lambda:self.show_page('Files')).pack(side='right')
        tk.Label(downloads_tab,text='Stop task pauses the current download. Queued and partial files are saved for resume after restarting the app.\nServers without a usable resume validator restart the partial transfer. Files appear in Files only after completion.',bg=BG,fg=MUTED,font=('DejaVu Sans',9),justify='left',wraplength=930).pack(anchor='w',pady=(0,10))
        self.person_analysis_table=self.tree(analyses_tab,[('label','Analysis',220),('sample','VCF sample',130),('build','Build',100),('date','Created',180),('path','Location',300)])
        self.analysis_entries={}
        for a in d['analyses']:
            self.analysis_entries[a['id']]=a;self.person_analysis_table.insert('','end',iid=a['id'],values=(a['label'],a['sample'],a['build'],a['created_at'],a['path']))
        self.person_analysis_table.bind('<Double-1>',lambda _:self.open_selected_analysis())
        row=tk.Frame(analyses_tab,bg=BG);row.pack(fill='x',pady=12);ttk.Button(row,text='Open selected analysis',command=self.open_selected_analysis).pack(side='left');ttk.Button(row,text='Import VCF analysis',command=self.import_dialog).pack(side='left',padx=10)
        tk.Label(self.content,text=str(self.person.path),bg=BG,fg=MUTED,font=('DejaVu Sans',8),anchor='w').pack(fill='x',pady=(12,0))
    def person_file_selected(self,_=None):
        ids=self.file_table.selection()
        if not ids:return
        r=self.file_entries[ids[0]];sha=r.get('sha256');integrity=('SHA-256 verified against supplied checksum: ' if r.get('checksum_verified') else 'Computed SHA-256: ')+sha if sha else 'No checksum recorded for this linked file.'
        self.file_note.configure(text=str(self.person.file_path(r))+'\n'+integrity)
    def copy_person_path(self):
        ids=self.file_table.selection()
        if ids:self.clipboard_clear();self.clipboard_append(str(self.person.file_path(self.file_entries[ids[0]])))
    def open_person_folder(self):
        if self.person:subprocess.Popen(['xdg-open',str(self.person.path)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    def add_person_files(self,copy=False):
        if self.busy or not self.person:return
        paths=filedialog.askopenfilenames(title='Copy files into this person project' if copy else 'Link existing files to this person',filetypes=[('Sequencing and analysis files','*.fastq *.fastq.gz *.fq *.fq.gz *.bam *.bai *.cram *.crai *.vcf *.vcf.gz *.tbi *.csi *.bed *.tsv *.csv *.pdf *.html *.md'),('All files','*')])
        if not paths:return
        owner=self.person
        def run():
            for path in paths:
                if self.cancel.is_set():raise Cancelled('File import stopped.')
                if copy:owner.copy_file(path,progress=lambda x:self.events.put(('progress',x)),cancel=self.cancel)
                else:owner.link_file(path)
            return len(paths)
        def done(n):self.show_page('Files');self.status.configure(text=f'{n} files '+('copied into' if copy else 'linked to')+' this person project.')
        self.start_task(run,done,'Adding files to person project…')
    def download_dialog(self):
        if self.busy or not self.person:return
        win=tk.Toplevel(self);win.title('Download into person project');win.geometry('680x485');win.configure(bg=BG);win.transient(self);win.grab_set()
        tk.Label(win,text='Download into '+self.person.data['name'],bg=BG,fg=INK,font=('DejaVu Sans',16,'bold'),wraplength=620).pack(anchor='w',padx=25,pady=(22,12))
        tk.Label(win,text='Direct HTTP / HTTPS file links — one link per line',bg=BG,fg=MUTED).pack(anchor='w',padx=25)
        links=ScrolledText(win,height=7,font=('DejaVu Sans',10),wrap='word');links.pack(fill='x',padx=25,pady=9)
        form=tk.Frame(win,bg=BG);form.pack(fill='x',padx=25,pady=8);filename=tk.StringVar();checksum=tk.StringVar()
        for row,(label,var) in enumerate([('Filename (optional, one link)',filename),('SHA-256 (optional, one link)',checksum)]):
            tk.Label(form,text=label,bg=BG,fg=MUTED,font=('DejaVu Sans',9)).grid(row=row,column=0,sticky='w',pady=7);ttk.Entry(form,textvariable=var,width=38).grid(row=row,column=1,padx=12,sticky='ew');form.columnconfigure(1,weight=1)
        tk.Label(win,text='Paired FASTQs can be queued together. Filenames normally come from each URL.\nDownloads are stored by type inside this person’s folder; existing files are never overwritten.',bg=BG,fg=MUTED,font=('DejaVu Sans',9),justify='left',wraplength=620).pack(anchor='w',padx=25,pady=10)
        def start():
            urls=[line.strip() for line in links.get('1.0','end').splitlines() if line.strip()]
            if not urls:messagebox.showerror('Download link needed','Paste at least one direct download link.',parent=win);return
            if len(urls)>1 and (filename.get().strip() or checksum.get().strip()):messagebox.showerror('One link required','Custom filename/checksum applies to one link. Leave those fields blank to queue multiple links.',parent=win);return
            try:jobs=[self.person.queue_download(url,filename.get().strip() or None,checksum.get()) for url in urls]
            except Exception as e:messagebox.showerror('Cannot queue download',str(e),parent=win);return
            win.destroy();self.start_download_jobs(jobs)
        ttk.Button(win,text='Download files',style='Accent.TButton',command=start).pack(anchor='e',padx=25,pady=12)
    def start_download_jobs(self,jobs):
        if self.busy or not self.person:return
        owner=self.person
        def run():
            result={'complete':0,'failed':0}
            for job in jobs:
                if self.cancel.is_set():raise Cancelled('Download queue paused; remaining jobs are saved.')
                try:owner.download(job['url'],job['filename'],job.get('expected_sha256',''),progress=lambda x:self.events.put(('progress',x)),cancel=self.cancel,job_id=job['id']);result['complete']+=1
                except Cancelled:raise
                except Exception:result['failed']+=1
            return result
        def done(result):self.show_page('Files');self.status.configure(text=f"Downloads: {result['complete']} complete, {result['failed']} failed. Open Downloads to retry failed jobs.")
        self.start_task(run,done,'Starting download queue…')
    def resume_selected_download(self):
        if self.busy:return
        ids=self.download_table.selection()
        if ids:self.start_download_jobs([self.download_entries[ids[0]]])
    def resume_pending_downloads(self):
        if self.busy or not self.person:return
        jobs=[j for j in self.person.downloads() if j['status']!='complete']
        if jobs:self.start_download_jobs(jobs)
        else:self.status.configure(text='No pending downloads for this person.')
    def use_person_file(self):
        if self.busy:return
        ids=self.file_table.selection()
        if not ids:return
        r=self.file_entries[ids[0]];path=self.person.file_path(r)
        if not path.is_file():messagebox.showerror('File unavailable','The linked file has moved or is missing.');return
        if path.name.lower().endswith(('.vcf','.vcf.gz')):self.import_dialog(path=path)
        elif path.suffix.lower()=='.bam':
            self.show_page('Run analysis');self.workflow_mode.set('Analysis-ready BAM');self.workflow_input.set(str(path))
        else:messagebox.showinfo('File stored in project','This file is organized in the project. Import a VCF or choose a prepared BAM for the analysis workflow. FASTQ alignment is not included yet.')
    def open_selected_analysis(self):
        if self.busy:return
        ids=self.person_analysis_table.selection()
        if ids:
            try:self.loaded(Project(self.analysis_entries[ids[0]]['path']),person=self.person)
            except Exception as e:messagebox.showerror('Analysis unavailable',str(e))
