import csv,gzip,json,sys,tempfile,threading,unittest
from pathlib import Path
from wes_lab.data import import_vcf,Project,classify_assertion,read_header,Cancelled
from wes_lab.pipeline import Pipeline

HEADER='''##fileformat=VCFv4.2
##reference=GRCh38
##INFO=<ID=POP_AF,Number=A,Type=Float,Description="Population allele frequency">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinVar assertion">
##INFO=<ID=ANN,Number=.,Type=String,Description="SnpEff annotation">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB
'''
def ann(allele,gene='PTEN',effect='missense_variant',impact='MODERATE',transcript='ENST00000371953.8'):
    return '|'.join([allele,effect,impact,gene,('ENSG00000161513' if gene=='FDXR' else 'ENSG00000171862'),'transcript',transcript,'protein_coding','1/9','c.1C>T','p.Ala1Val','','','','',''])

class ImportTests(unittest.TestCase):
    def setUp(self):self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()
    def vcf(self,records,name='input.vcf',header=HEADER):
        p=self.root/name;text=header+records
        if name.endswith('.gz'):
            with gzip.open(p,'wt') as f:f.write(text)
        else:p.write_text(text)
        return p
    def load(self,records,**kwargs):
        path=self.vcf(records);return import_vcf(path,self.root/'project',sample='A',build='GRCh38',**kwargs)
    def test_conflict_is_not_a_pathogenic_assertion(self):
        self.assertEqual(classify_assertion('Conflicting_classifications_of_pathogenicity'),'Conflicting')
        self.assertEqual(classify_assertion('Pathogenic/Likely_pathogenic'),'Pathogenic / likely pathogenic')
        p=self.load('1\t12\t.\tC\tT\t50\tPASS\tCLNSIG=Conflicting_classifications_of_pathogenicity;ANN='+ann('T')+'\tGT:DP:GQ:AD\t0/1:30:99:15,15\t0/0:20:60:20,0\n')
        self.assertEqual(p.summary()['clinical_assertions'],0)
        self.assertEqual(p.variants()[1][0]['review_group'],'Conflicting assertions')
    def test_multiallelic_sample_and_frequency_are_allele_specific(self):
        text='1\t12\t.\tC\tT,G\t50\tPASS\tPOP_AF=0.0001,0.2;CLNSIG=Pathogenic;ANN='+ann('T',impact='HIGH')+','+ann('G',gene='FDXR',effect='synonymous_variant',impact='LOW')+'\tGT:DP:GQ:AD\t0/2:20:99:10,0,10\t0/1:20:99:10,10,0\n'
        p=self.load(text);count,rows=p.variants();self.assertEqual(count,1);r=rows[0]
        self.assertEqual(r['alt'],'G');self.assertEqual(r['genotype'],'0/2');self.assertAlmostEqual(r['alt_fraction'],.5);self.assertAlmostEqual(r['population_af'],.2)
        self.assertEqual(r['effect'],'synonymous_variant');self.assertEqual(r['review_group'],'Support check')
        self.assertEqual(p.detail(r['id'])['annotations'][0]['gene'],'FDXR')
    def test_reference_genotypes_are_not_sample_variants(self):
        p=self.load('1\t10\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ\t0/0:20:99\t0/1:20:99\n')
        self.assertEqual(p.summary()['total'],0)
    def test_mane_preferred_over_alternative_frameshift(self):
        text='10\t12\t.\tC\tT\t50\tPASS\tANN='+ann('T',effect='frameshift_variant',impact='HIGH',transcript='FAKE.1')+','+ann('T',effect='5_prime_UTR_variant',impact='MODIFIER')+'\tGT:DP:GQ:AD\t0/1:30:99:15,15\t0/0:20:60:20,0\n'
        p=self.load(text);r=p.variants()[1][0];self.assertEqual(r['impact'],'MODIFIER');self.assertEqual(r['review_group'],'Transcript check');self.assertEqual(len(p.detail(r['id'])['annotations']),2)
    def test_generic_af_is_not_population_af(self):
        p=self.load('1\t12\t.\tC\tT\t50\tPASS\tAF=0.5;ANN='+ann('T')+'\tGT:DP:GQ:AD\t0/1:30:99:15,15\t0/0:20:60:20,0\n')
        self.assertIsNone(p.variants()[1][0]['population_af'])
    def test_gzip_and_required_multisample_selection(self):
        path=self.vcf('1\t12\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:99\t0/0:20:60\n','input.vcf.gz')
        self.assertEqual(read_header(path)['samples'],['A','B'])
        with self.assertRaises(ValueError):import_vcf(path,self.root/'new')
        self.assertEqual(import_vcf(path,self.root/'new',sample='B').summary()['total'],0)
    def test_failed_reimport_preserves_existing_project(self):
        p=self.load('1\t12\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:99\t0/0:20:60\n')
        bad=self.vcf('1\tBROKEN\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\t0/0\n','bad.vcf')
        with self.assertRaises(ValueError):import_vcf(bad,p.path,sample='A')
        self.assertEqual(Project(p.path).summary()['total'],1);self.assertFalse(list(p.path.glob('.import-*')))
    def test_filter_and_html_export_handle_untrusted_text(self):
        path=self.vcf('1\t12\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:99\t0/0:20:60\n',header=HEADER.replace('\tA\tB\n','\t<script>\tB\n'))
        p=import_vcf(path,self.root/'project',sample='<script>');self.assertEqual(p.variants({'query':"' OR 1=1 --"})[0],0)
        out=self.root/'report.html';p.export_report(out);self.assertNotIn('<script>',out.read_text());self.assertIn('&lt;script&gt;',out.read_text())
    def test_bad_coverage_does_not_replace_valid_coverage(self):
        p=self.load('1\t12\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:99\t0/0:20:60\n')
        cov=self.root/'coverage.tsv';cov.write_text('gene\tmean_depth\tpercent_ge20\nPTEN\t40\t95\n');p.import_coverage(cov)
        cov.write_text('gene\tmean_depth\tpercent_ge20\nPTEN\t40\t150\n')
        with self.assertRaises(ValueError):p.import_coverage(cov)
        self.assertEqual(p.summary()['coverage_complete'],1)
    def test_vep_annotation_matches_alt_number(self):
        header=HEADER.replace('##INFO=<ID=ANN,Number=.,Type=String,Description="SnpEff annotation">','##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|HGVSc|HGVSp|ALLELE_NUM|MANE_SELECT">')
        text='1\t12\t.\tC\tT,G\t50\tPASS\tCSQ=T|stop_gained|HIGH|PTEN|GENE1|TR1|c.1|p.1|1|NM_1,G|synonymous_variant|LOW|FDXR|GENE2|TR2|c.2|p.2|2|NM_2\tGT:DP:GQ\t0/2:30:99\t0/1:30:99\n'
        path=self.vcf(text,header=header);p=import_vcf(path,self.root/'project',sample='A');r=p.variants()[1][0];self.assertEqual(r['genes'],'FDXR');self.assertEqual(r['impact'],'LOW');self.assertEqual(r['mane'],1)
    def test_cancel_import_publishes_no_partial_project(self):
        event=threading.Event();event.set();path=self.vcf('1\t12\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\t0/0\n')
        with self.assertRaises(Cancelled):import_vcf(path,self.root/'project',sample='A',cancel=event)
        self.assertFalse((self.root/'project/variants.sqlite').exists())
    def test_spanning_deletion_placeholders_are_counted_but_not_presented_as_alleles(self):
        p=self.load('1\t12\t.\tC\t*\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:99\t0/0:20:60\n')
        self.assertEqual(p.summary()['total'],0)
        self.assertEqual(p.meta['counts']['input_records'],1)
        self.assertEqual(p.meta['counts']['placeholder_alleles_skipped'],1)
    def test_pipeline_stop_terminates_its_process_group(self):
        event=threading.Event();pipeline=Pipeline(self.root/'input.bam',self.root/'workflow',cancel=event);timer=threading.Timer(.25,event.set);timer.start()
        try:
            with self.assertRaises(Cancelled):pipeline.execute([sys.executable,'-c','import time; time.sleep(10)'])
        finally:timer.cancel()
        self.assertIsNone(pipeline.process);self.assertTrue((pipeline.out/'workflow.log').exists())

if __name__=='__main__':unittest.main()
