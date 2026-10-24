"""One person per project, with managed files, linked analyses and resumable downloads."""
from __future__ import annotations
import contextlib,hashlib,json,os,re,threading,time,uuid
from pathlib import Path
from urllib.parse import urlsplit,urlunsplit,unquote
from urllib.request import Request,build_opener,HTTPRedirectHandler
from urllib.error import HTTPError
from .data import Cancelled

CHUNK=1024*1024
FOLDERS=('files/fastq','files/bam','files/vcf','files/other','analyses','reports','downloads')
def now():return time.strftime('%Y-%m-%d %H:%M:%S')
def clean_name(name):
    if not name or name in ('.','..') or '/' in name or '\\' in name or any(ord(c)<32 for c in name):raise ValueError('Use a filename without folder paths or control characters.')
    if len(name.encode())>230:raise ValueError('Filename is too long.')
    return name

def file_kind(name):
    n=name.lower()
    if n.endswith(('.fastq','.fastq.gz','.fq','.fq.gz')):return 'FASTQ','files/fastq'
    if n.endswith(('.bam','.bai','.cram','.crai')):return ('Alignment index' if n.endswith(('.bai','.crai')) else 'BAM / CRAM'),'files/bam'
    if n.endswith(('.vcf','.vcf.gz','.bcf','.tbi','.csi')):return ('Variant index' if n.endswith(('.tbi','.csi')) else 'VCF / BCF'),'files/vcf'
    if n.endswith(('.html','.pdf','.md')):return 'Report','reports'
    return 'Other','files/other'

def public_url(url):
    p=urlsplit(url);return urlunsplit((p.scheme,p.netloc.split('@')[-1],p.path,'',''))

def validate_url(url):
    p=urlsplit(url)
    if p.scheme not in ('http','https') or not p.hostname:raise ValueError('Paste a direct HTTP or HTTPS download link.')
    if p.username or p.password:raise ValueError('Use an HTTPS signed link instead of embedding a username/password in the URL.')
    if any(ord(c)<32 for c in url):raise ValueError('Download URL contains invalid characters.')
    return url

class SafeRedirect(HTTPRedirectHandler):
    def redirect_request(self,req,fp,code,msg,headers,newurl):
        validate_url(newurl)
        if urlsplit(req.full_url).scheme=='https' and urlsplit(newurl).scheme!='https':raise ValueError('The server tried to redirect an HTTPS download to an insecure connection.')
        return super().redirect_request(req,fp,code,msg,headers,newurl)

def atomic_json(path,data):
    path=Path(path);temporary=path.with_name('.'+path.name+'.'+uuid.uuid4().hex+'.tmp')
    try:
        with temporary.open('x',encoding='utf-8') as f:
            os.chmod(temporary,0o600);json.dump(data,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
        os.replace(temporary,path)
    finally:temporary.unlink(missing_ok=True)

@contextlib.contextmanager
def lock_file(path):
    import fcntl
    with Path(path).open('a+') as lock:
        os.chmod(path,0o600)
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError('This project operation is already running. Try again when it finishes.') from None
        yield

class PersonProject:
    def __init__(self,path):
        self.path=Path(path).resolve();self.manifest=self.path/'person.json'
        self.data=self.read()
        if self.data.get('schema_version')!=1 or not self.data.get('id'):raise ValueError('Not a supported person project.')
    @classmethod
    def create(cls,parent,name,notes=''):
        name=name.strip()
        if not name or len(name)>150 or any(ord(c)<32 for c in name):raise ValueError('Enter a person name or identifier (up to 150 characters).')
        parent=Path(parent).resolve();parent.mkdir(parents=True,exist_ok=True)
        slug=re.sub(r'[^A-Za-z0-9._-]+','_',name).strip('._')[:55] or 'person'
        path=parent/(slug+'_'+uuid.uuid4().hex[:8]);path.mkdir(mode=0o700)
        for folder in FOLDERS:(path/folder).mkdir(parents=True,exist_ok=True,mode=0o700)
        atomic_json(path/'person.json',dict(schema_version=1,id=uuid.uuid4().hex,name=name,notes=notes,created_at=now(),updated_at=now(),files=[],analyses=[],active_analysis=None))
        return cls(path)
    def read(self):
        try:return json.loads(self.manifest.read_text())
        except (OSError,ValueError) as e:raise ValueError(f'Cannot open person project: {self.path.name}') from e
    def change(self,fn):
        with lock_file(self.path/'.person.lock'):
            data=self.read();result=fn(data);data['updated_at']=now();atomic_json(self.manifest,data);self.data=data;return result
    def update(self,name,notes):
        name=name.strip()
        if not name or len(name)>150 or any(ord(c)<32 for c in name):raise ValueError('Enter a valid name or person identifier.')
        def mutate(d):d.update(name=name,notes=notes)
        self.change(mutate)
    def inside(self,relative):
        path=(self.path/relative).resolve()
        if not path.is_relative_to(self.path) or path==self.path:raise ValueError('Project paths must stay inside this person’s folder.')
        return path
    def destination(self,name):
        name=clean_name(name);_,folder=file_kind(name);return self.inside(Path(folder)/name)
    def file_path(self,entry):
        return self.inside(entry['path']) if entry['storage']=='managed' else Path(entry['path']).resolve()
    def register_file(self,path,kind=None,source='',sha256=None,checksum_verified=False):
        path=Path(path).resolve()
        if not path.is_file():raise ValueError('Select an existing file.')
        managed=path.is_relative_to(self.path);kind=kind or file_kind(path.name)[0]
        entry=dict(id=uuid.uuid4().hex,name=path.name,path=str(path.relative_to(self.path)) if managed else str(path),storage='managed' if managed else 'linked',kind=kind,size=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns,added_at=now(),source=source,sha256=sha256,checksum_verified=checksum_verified)
        def mutate(d):
            for old in d['files']:
                if self.file_path(old)==path:return old
            d['files'].append(entry);return entry
        return self.change(mutate)
    def link_file(self,path):return self.register_file(path,source='Linked local file; not copied')
    def copy_file(self,source,progress=None,cancel=None):
        source=Path(source).resolve()
        if not source.is_file():raise ValueError('Source file is missing.')
        before=source.stat()
        dest=self.destination(source.name)
        if source==dest:return self.register_file(dest)
        if dest.exists():raise FileExistsError(f'{dest.name} already exists in this project. No file was overwritten.')
        partial=self.inside('downloads/copy_'+uuid.uuid4().hex+'.part');digest=hashlib.sha256();done=0
        try:
            with source.open('rb') as src,partial.open('xb') as out:
                os.chmod(partial,0o600)
                for block in iter(lambda:src.read(CHUNK),b''):
                    if cancel and cancel.is_set():raise Cancelled('Copy stopped; source file unchanged.')
                    out.write(block);digest.update(block);done+=len(block)
                    if progress:progress(f'Copying {source.name}: {done:,} / {source.stat().st_size:,} bytes')
                out.flush();os.fsync(out.fileno())
            after=source.stat()
            if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):raise ValueError('Source file changed during copy. No completed copy was published.')
            os.link(partial,dest) # Atomic, fails if another operation already created dest.
            return self.register_file(dest,source='Copied local file',sha256=digest.hexdigest())
        finally:partial.unlink(missing_ok=True)
    def files(self):
        rows=[]
        for entry in self.read()['files']:
            r=dict(entry);p=self.file_path(r);r['exists']=p.is_file();r['changed']=r['exists'] and (p.stat().st_size!=r['size'] or p.stat().st_mtime_ns!=r['mtime_ns']);rows.append(r)
        return rows
    def new_analysis_directory(self,name='analysis'):
        slug=re.sub(r'[^A-Za-z0-9._-]+','_',name).strip('._')[:60] or 'analysis'
        path=self.inside(Path('analyses')/(slug+'_'+time.strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:6]));path.mkdir();return path
    def register_analysis(self,path,sample,build,label=None):
        path=Path(path).resolve()
        if not (path/'variants.sqlite').is_file():raise ValueError('Analysis database is missing.')
        def mutate(d):
            for a in d['analyses']:
                if Path(a['path']).resolve()==path:d['active_analysis']=a['id'];return a
            a=dict(id=uuid.uuid4().hex,path=str(path),sample=sample,build=build,label=label or sample,created_at=now());d['analyses'].append(a);d['active_analysis']=a['id'];return a
        return self.change(mutate)
    def downloads(self):
        rows=[]
        for p in self.inside('downloads').glob('*.json'):
            try:
                row=json.loads(p.read_text());row['source_display']=public_url(row['url']);partial=self.inside('downloads/'+row['id']+'.part');row['bytes']=partial.stat().st_size if partial.exists() else row.get('bytes',0);rows.append(row)
            except (OSError,ValueError,KeyError):continue
        return sorted(rows,key=lambda r:r.get('updated_at',''),reverse=True)
    def queue_download(self,url,filename=None,expected_sha256=''):
        url=validate_url(url.strip());filename=clean_name(filename or unquote(Path(urlsplit(url).path).name) or 'download.bin')
        expected_sha256=expected_sha256.strip().lower()
        if expected_sha256 and not re.fullmatch('[0-9a-f]{64}',expected_sha256):raise ValueError('SHA-256 must be 64 hexadecimal characters.')
        ident=hashlib.sha256((url+'\0'+filename).encode()).hexdigest()[:24];path=self.inside('downloads/'+ident+'.json')
        with lock_file(self.inside('downloads/'+ident+'.lock')):
            if path.exists():
                state=json.loads(path.read_text())
                if expected_sha256:state['expected_sha256']=expected_sha256;atomic_json(path,state)
            else:
                if self.destination(filename).exists():raise FileExistsError(filename+' already exists in this project.')
                state=dict(id=ident,url=url,filename=filename,expected_sha256=expected_sha256,status='queued',bytes=0,total=None,created_at=now(),updated_at=now());atomic_json(path,state)
        return state
    def download(self,url,filename=None,expected_sha256='',progress=None,cancel=None,job_id=None):
        url=validate_url(url.strip());filename=clean_name(filename or unquote(Path(urlsplit(url).path).name) or 'download.bin');expected_sha256=expected_sha256.strip().lower()
        if expected_sha256 and not re.fullmatch('[0-9a-f]{64}',expected_sha256):raise ValueError('SHA-256 must be 64 hexadecimal characters, or left blank.')
        ident=hashlib.sha256((url+'\0'+filename).encode()).hexdigest()[:24]
        if job_id and job_id!=ident:raise ValueError('Download identity does not match the saved job.')
        dest=self.destination(filename);partial=self.inside('downloads/'+ident+'.part');state_path=self.inside('downloads/'+ident+'.json')
        with lock_file(self.inside('downloads/'+ident+'.lock')):
            state=json.loads(state_path.read_text()) if state_path.exists() else dict(id=ident,url=url,filename=filename,status='queued',bytes=0,total=None,created_at=now())
            if dest.exists():
                if state.get('status')=='complete' and state.get('sha256'):
                    existing=next((r for r in self.files() if self.file_path(r)==dest),None)
                    if existing and not existing['changed']:
                        if expected_sha256 and expected_sha256!=state['sha256']:raise ValueError('Existing file does not match the requested SHA-256.')
                        return existing
                state.update(status='failed',error='Destination filename already exists; choose another filename.',updated_at=now());atomic_json(state_path,state)
                raise FileExistsError(f'{filename} already exists. Choose a different filename; no file was overwritten.')
            state['expected_sha256']=expected_sha256 or state.get('expected_sha256','')
            def save(status):state.update(status=status,bytes=partial.stat().st_size if partial.exists() else state.get('bytes',0),updated_at=now());atomic_json(state_path,state)
            def check_cancel():
                if cancel and cancel.is_set():raise Cancelled('Download paused. Select it in Files and choose Resume download.')
            opener=build_opener(SafeRedirect());response=None
            try:
                check_cancel();offset=partial.stat().st_size if partial.exists() else 0
                validator=state.get('etag')
                if validator and validator.startswith('W/'):validator=None
                validator=validator or state.get('last_modified')
                # Without a validator, do not append bytes from a potentially changed object.
                resume=bool(offset and validator)
                headers={'User-Agent':'WES-Lab/0.2','Accept-Encoding':'identity'}
                if resume:headers.update(Range=f'bytes={offset}-',**{'If-Range':validator})
                save('connecting')
                try:response=opener.open(Request(url,headers=headers),timeout=30)
                except HTTPError as e:
                    if e.code!=416 or not resume:raise
                    response=opener.open(Request(url,headers={'User-Agent':'WES-Lab/0.2','Accept-Encoding':'identity'}),timeout=30);resume=False
                with response:
                    validate_url(response.geturl())
                    if response.headers.get('Content-Encoding','identity').lower()!='identity':raise ValueError('Server returned encoded content despite an identity request; cannot safely resume.')
                    if 'text/html' in response.headers.get('Content-Type','').lower():raise ValueError('This link returned a web page. Paste the direct file-download URL instead.')
                    append=resume and response.status==206
                    content_range=response.headers.get('Content-Range','')
                    if response.status==206:
                        match=re.fullmatch(r'bytes (\d+)-(\d+)/(\d+|\*)',content_range)
                        if not append or not match or int(match[1])!=offset:raise ValueError('Server returned an invalid resumed range; partial file retained.')
                        if state.get('etag') and response.headers.get('ETag') and state['etag']!=response.headers['ETag']:raise ValueError('Remote file changed during resume. Use a new filename to restart safely.')
                        total=int(match[3]) if match[3]!='*' else None
                    elif response.status==200:
                        offset=0;append=False;total=int(response.headers['Content-Length']) if response.headers.get('Content-Length') else None
                    else:raise ValueError(f'Unexpected HTTP status {response.status}.')
                    state.update(etag=response.headers.get('ETag'),last_modified=response.headers.get('Last-Modified'),total=total)
                    with partial.open('ab' if append else 'wb') as out:
                        os.chmod(partial,0o600);done=offset;last_save=0
                        for block in iter(lambda:response.read(CHUNK),b''):
                            check_cancel();out.write(block);done+=len(block);state['bytes']=done
                            if progress:progress(f'Downloading {filename}: {done:,}'+(f' / {total:,} bytes ({100*done/total:.1f}%)' if total else ' bytes'))
                            if time.monotonic()-last_save>1:out.flush();save('downloading');last_save=time.monotonic()
                        out.flush();os.fsync(out.fileno())
                    if total is not None and done!=total:raise ValueError(f'Download is incomplete ({done:,} of {total:,} bytes). Resume to continue.')
                check_cancel();save('verifying')
                if progress:progress(f'Verifying SHA-256 for {filename}…')
                digest=hashlib.sha256()
                with partial.open('rb') as src:
                    for block in iter(lambda:src.read(CHUNK),b''):check_cancel();digest.update(block)
                sha=digest.hexdigest()
                if state['expected_sha256'] and sha!=state['expected_sha256']:raise ValueError('SHA-256 does not match the supplied checksum. The file was not published.')
                check_cancel();os.link(partial,dest)
                result=self.register_file(dest,source=public_url(url),sha256=sha,checksum_verified=bool(state['expected_sha256']))
                state.update(sha256=sha,bytes=dest.stat().st_size);save('complete');partial.unlink(missing_ok=True);return result
            except Cancelled:save('paused');raise
            except BaseException as e:state['error']=str(e);save('failed');raise

def list_people(parent):
    rows=[]
    for manifest in Path(parent).glob('*/person.json'):
        try:
            person=PersonProject(manifest.parent);rows.append(person)
        except ValueError:continue
    return sorted(rows,key=lambda p:p.data.get('updated_at',''),reverse=True)

def find_analysis_owner(parent,analysis):
    analysis=Path(analysis).resolve()
    for person in list_people(parent):
        if any(Path(a['path']).resolve()==analysis for a in person.data['analyses']):return person
    return None
