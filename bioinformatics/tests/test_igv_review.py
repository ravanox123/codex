import unittest
import xml.etree.ElementTree as ET
from wes_lab.igv_review import locus,prepare_session

class IGVReviewTests(unittest.TestCase):
    def test_reference_alias_and_bounds(self):
        self.assertEqual(locus('17',74863914),'CM000679.2:74863814-74864014')
        self.assertEqual(locus('chr2',1),'CM000664.2:1-101')
        with self.assertRaises(ValueError):locus('invalid',12)
        with self.assertRaises(ValueError):locus('2',0)
    def test_session_local_resources_and_indexes(self):
        from pathlib import Path
        path=prepare_session();root=ET.parse(path).getroot()
        self.assertTrue(Path(root.attrib['genome']+'.fai').exists())
        resources=root.findall('Resources/Resource')
        self.assertEqual(len(resources),4)
        self.assertTrue(all(Path(r.attrib['path']).is_file() for r in resources))
        self.assertEqual(root.attrib['locus'],locus('17',74863914))

if __name__=='__main__':unittest.main()
