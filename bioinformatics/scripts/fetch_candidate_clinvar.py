"""Cache public ClinVar records for known candidate IDs, without uploading sample files."""
import json, time, urllib.request
from pathlib import Path
from html.parser import HTMLParser
IDS = dict(USH2A=48427,EXOC6B=1601386,DST=809956,ENPP1=13589,DLX6=1385415,STXBP1=365052,NOTCH1=477890,POLR3A=368925,PTEN=440216,SYCE1=2672423,CTR9=1091039,VDR=1252227,BRCA2=51231,MYO1E=983002,ZFHX3=808107,POLR2A=1686168,SETBP1=326725,KATNAL2=3060277,FDXR=520994,WNT10A=464181,VWF=296,DNAH14=4880149,CNN2_1037757=1696852,CNN2_1037767=1696851,CDKN2B=812639,TRMT5=372246,PSPH=403351)
OUT=Path(__file__).resolve().parents[1]/'results/wgs/candidate_review'
class Text(HTMLParser):
 def __init__(self):super().__init__();self.parts=[];self.skip=0
 def handle_starttag(self,tag,attrs):
  if tag in ('script','style'):self.skip+=1
 def handle_endtag(self,tag):
  if tag in ('script','style'):self.skip=max(0,self.skip-1)
 def handle_data(self,data):
  if not self.skip and data.strip():self.parts.append(data.strip())
for gene,vid in IDS.items():
 p=OUT/f'clinvar_{vid}.html';url=f'https://www.ncbi.nlm.nih.gov/clinvar/variation/{vid}/'
 try:
  if not p.exists():
   req=urllib.request.Request(url,headers={'User-Agent':'LocalResearchVariantReview/1.0'})
   p.write_bytes(urllib.request.urlopen(req,timeout=45).read());time.sleep(.4)
  parser=Text();parser.feed(p.read_text());text='\n'.join(parser.parts);(OUT/f'clinvar_{vid}.txt').write_text(text)
  a=text.find('Germline\nClassification');b=text.find('Somatic',a)
  print(gene,vid,text[a:b].replace('\n',' ')[:650],flush=True)
 except Exception as e:print(gene,vid,'FETCH_FAILED',str(e),flush=True)
