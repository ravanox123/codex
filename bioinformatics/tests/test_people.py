import hashlib,json,tempfile,threading,unittest
from pathlib import Path
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from urllib.error import HTTPError
from wes_lab.people import PersonProject,list_people,find_analysis_owner,CHUNK,public_url
from wes_lab.data import Cancelled

class PersonTests(unittest.TestCase):
    def setUp(self):self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name);self.person=PersonProject.create(self.root/'people','Person A','Notes')
    def tearDown(self):self.tmp.cleanup()
    def test_separate_persons_with_same_display_name_are_isolated(self):
        other=PersonProject.create(self.root/'people','Person A');self.assertNotEqual(other.path,self.person.path)
        self.assertTrue((self.person.path/'files/fastq').is_dir());self.assertEqual(len(list_people(self.root/'people')),2)
        self.person.update('New label','Updated notes');self.assertEqual(PersonProject(self.person.path).data['notes'],'Updated notes')
    def test_link_and_copy_keep_originals_and_detect_missing_or_changed_files(self):
        source=self.root/'reads.fastq';source.write_bytes(b'ACGT')
        linked=self.person.link_file(source);self.assertEqual(linked['storage'],'linked');self.assertEqual(len(self.person.files()),1)
        self.person.link_file(source);self.assertEqual(len(self.person.files()),1)
        copied=self.person.copy_file(source);self.assertEqual(copied['storage'],'managed');self.assertEqual(self.person.file_path(copied).read_bytes(),b'ACGT')
        source.write_bytes(b'CHANGED');self.assertTrue(next(x for x in self.person.files() if x['id']==linked['id'])['changed'])
        source.unlink();self.assertFalse(next(x for x in self.person.files() if x['id']==linked['id'])['exists']);self.assertTrue(self.person.file_path(copied).exists())
    def test_copy_cancellation_publishes_no_partial_file(self):
        source=self.root/'sample.bam';source.write_bytes(b'A'*(CHUNK*2));cancel=threading.Event()
        with self.assertRaises(Cancelled):self.person.copy_file(source,progress=lambda _:cancel.set(),cancel=cancel)
        self.assertFalse(self.person.destination(source.name).exists());self.assertFalse(list((self.person.path/'downloads').glob('*.part')))
        self.assertEqual(source.stat().st_size,CHUNK*2)
    def test_filename_and_symlink_escape_are_rejected(self):
        for name in ('../sample.bam','/tmp/sample.bam','..','bad\\name.fastq'):
            with self.assertRaises(ValueError):self.person.destination(name)
        (self.person.path/'files/bam').rmdir();(self.person.path/'files/bam').symlink_to(self.root)
        with self.assertRaises(ValueError):self.person.destination('sample.bam')
    def test_analysis_ownership_survives_reopening(self):
        path=self.person.new_analysis_directory('exome');(path/'variants.sqlite').touch();self.person.register_analysis(path,'sample-A','GRCh38')
        self.person.register_analysis(path,'sample-A','GRCh38');self.assertEqual(len(self.person.read()['analyses']),1)
        self.assertEqual(find_analysis_owner(self.root/'people',path).data['id'],self.person.data['id'])
    def test_urls_are_validated_and_display_does_not_reveal_signed_queries(self):
        for url in ('file:///etc/passwd','ftp://host/a.bam','https://user:password@host/a.bam'):
            with self.assertRaises(ValueError):self.person.queue_download(url)
        job=self.person.queue_download('https://example.org/reads.fastq?token=secret')
        row=self.person.downloads()[0];self.assertNotIn('secret',row['source_display']);self.assertEqual(row['status'],'queued')
        self.assertEqual(public_url(job['url']),'https://example.org/reads.fastq')

PAYLOAD=(b'@read\nACGTACGT\n+\nIIIIIIII\n'*(CHUNK//24*3))
class Handler(BaseHTTPRequestHandler):
    calls=[]
    def log_message(self,*args):pass
    def do_GET(self):
        type(self).calls.append((self.path,self.headers.get('Range'),self.headers.get('If-Range')))
        if self.path=='/html':
            self.send_response(200);self.send_header('Content-Type','text/html');self.end_headers();self.wfile.write(b'<html>Login</html>');return
        if self.path=='/redirect':
            self.send_response(302);self.send_header('Location','file:///etc/passwd');self.end_headers();return
        resume=self.path not in ('/noresume','/truncated') and self.headers.get('Range')
        start=int(self.headers['Range'].split('=')[1].split('-')[0]) if resume else 0
        if start>=len(PAYLOAD):self.send_response(416);self.end_headers();return
        self.send_response(206 if resume else 200);self.send_header('Content-Type','application/octet-stream')
        if self.path!='/noresume':self.send_header('ETag','"v1"')
        self.send_header('Content-Length',str(len(PAYLOAD)-start+(20 if self.path=='/truncated' else 0)))
        if resume:self.send_header('Content-Range',f'bytes {start+(1 if self.path=="/bad_range" else 0)}-{len(PAYLOAD)-1}/{len(PAYLOAD)}')
        self.end_headers()
        try:self.wfile.write(PAYLOAD[start:]);self.wfile.flush()
        except (BrokenPipeError,ConnectionResetError):pass
        self.close_connection=True

class DownloadTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server=ThreadingHTTPServer(('127.0.0.1',0),Handler);cls.thread=threading.Thread(target=cls.server.serve_forever,daemon=True);cls.thread.start();cls.base=f'http://127.0.0.1:{cls.server.server_port}'
    @classmethod
    def tearDownClass(cls):cls.server.shutdown();cls.server.server_close();cls.thread.join()
    def setUp(self):self.tmp=tempfile.TemporaryDirectory();self.person=PersonProject.create(Path(self.tmp.name)/'people','Local download test');Handler.calls=[]
    def tearDown(self):self.tmp.cleanup()
    def pause(self,endpoint='/reads.fastq',filename='reads.fastq'):
        event=threading.Event();job=self.person.queue_download(self.base+endpoint,filename)
        with self.assertRaises(Cancelled):self.person.download(job['url'],job['filename'],progress=lambda _:event.set(),cancel=event)
        self.assertEqual(self.person.downloads()[0]['status'],'paused');self.assertFalse(self.person.destination(filename).exists());return job
    def test_download_pause_resume_preserves_exact_content_and_hash(self):
        job=self.pause();before=self.person.downloads()[0]['bytes'];self.assertGreater(before,0)
        digest=hashlib.sha256(PAYLOAD).hexdigest();entry=self.person.download(job['url'],job['filename'],digest)
        self.assertEqual(self.person.file_path(entry).read_bytes(),PAYLOAD);self.assertTrue(entry['checksum_verified'])
        self.assertIn(('/reads.fastq',f'bytes={before}-','"v1"'),Handler.calls)
        self.assertEqual(self.person.downloads()[0]['status'],'complete');self.assertFalse(list((self.person.path/'downloads').glob('*.part')))
    def test_server_without_validator_restarts_without_duplicating_bytes(self):
        job=self.pause('/noresume');entry=self.person.download(job['url'],job['filename']);self.assertEqual(self.person.file_path(entry).read_bytes(),PAYLOAD)
        self.assertTrue(all(call[1] is None for call in Handler.calls))
    def test_checksum_mismatch_never_publishes_file(self):
        with self.assertRaises(ValueError):self.person.download(self.base+'/reads.fastq','reads.fastq','0'*64)
        self.assertFalse(self.person.destination('reads.fastq').exists());self.assertEqual(self.person.downloads()[0]['status'],'failed')
    def test_existing_files_are_never_overwritten(self):
        dest=self.person.destination('sample.bam');dest.write_bytes(b'original')
        with self.assertRaises(FileExistsError):self.person.download(self.base+'/reads.fastq','sample.bam')
        self.assertEqual(dest.read_bytes(),b'original')
    def test_invalid_resumed_range_does_not_append_bytes(self):
        job=self.pause('/bad_range');partial=self.person.inside('downloads/'+job['id']+'.part');size=partial.stat().st_size
        with self.assertRaises(ValueError):self.person.download(job['url'],job['filename'])
        self.assertEqual(partial.stat().st_size,size);self.assertFalse(self.person.destination(job['filename']).exists())
    def test_html_incomplete_and_non_http_redirects_are_not_published(self):
        for endpoint in ('/html','/truncated','/redirect'):
            with self.assertRaises(Exception):self.person.download(self.base+endpoint,endpoint[1:]+'.fastq')
            self.assertFalse(self.person.destination(endpoint[1:]+'.fastq').exists())
    def test_queued_jobs_survive_reopening(self):
        self.person.queue_download(self.base+'/reads.fastq','R1.fastq');self.person.queue_download(self.base+'/reads.fastq','R2.fastq')
        reopened=PersonProject(self.person.path);self.assertEqual(len(reopened.downloads()),2)
        for job in reopened.downloads():reopened.download(job['url'],job['filename'])
        self.assertEqual(len(reopened.files()),2)

if __name__=='__main__':unittest.main()
