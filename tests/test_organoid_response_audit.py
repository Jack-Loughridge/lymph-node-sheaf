import csv
import math
from pathlib import Path
import sys
import tempfile
import unittest
from zipfile import ZipFile
from xml.sax.saxutils import escape

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src/organoid_audit'))
import audit_organoid_response_genes as a


def workbook(path, sheetrows):
    ns = a.NS['m']
    with ZipFile(path,'w') as z:
        z.writestr('xl/workbook.xml',f'<workbook xmlns="{ns}" xmlns:r="{a.NS["r"]}"><sheets>'+''.join(
            f'<sheet name="{name}" sheetId="{i}" r:id="r{i}"/>' for i,name in enumerate(sheetrows,1))+'</sheets></workbook>')
        z.writestr('xl/_rels/workbook.xml.rels','<Relationships>'+''.join(
            f'<Relationship Id="r{i}" Target="worksheets/sheet{i}.xml"/>' for i in range(1,len(sheetrows)+1))+'</Relationships>')
        for i,rows in enumerate(sheetrows.values(),1):
            xml = f'<worksheet xmlns="{ns}"><sheetData>'
            for ri,row in enumerate(rows,1):
                xml += f'<row r="{ri}">'
                for ci,value in enumerate(row):
                    ref=f'{chr(65+ci)}{ri}'
                    if isinstance(value,str):
                        xml+=f'<c r="{ref}" t="inlineStr"><is><t>{escape(value)}</t></is></c>'
                    elif value is not None:
                        xml+=f'<c r="{ref}"><v>{value}</v></c>'
                xml+='</row>'
            z.writestr(f'xl/worksheets/sheet{i}.xml',xml+'</sheetData></worksheet>')


class AuditTests(unittest.TestCase):
    def test_aliases(self):
        for label in ['TNFA','TNF-alpha','TNFALPHA','tnf']:
            self.assertEqual(a.factor_name(label),'TNF')
        self.assertEqual(a.factor_name('ISlR'),'ISLR')

    def test_strict_thresholds(self):
        for p,fc in [(0.05,2),(.01,1),(.01,-1),(-.1,2),(.01,math.nan),(math.nan,2)]:
            self.assertFalse(a.significant(p,fc))
        self.assertTrue(a.significant(0,-1.01))
        self.assertTrue(a.significant(.049,1.01))

    def test_duplicate_sheets_signed_genes_and_skips(self):
        header=['gene','baseMean','log2FoldChange','lfcSE','pvalue','padj','pert']
        duplicate=[header,['GENE1',1,2,1,.001,.01,'TNFA']]
        rows=[header,['GENE2',1,-2,1,.001,.01,'TNF-alpha'],
              ['GENE3',1,1,1,.001,.01,'TNFA'],
              ['GENE4',1,2,1,.001,.05,'TNFA'],
              ['GENE5',1,2,1,.001,None,'TNFA'],
              ['GENE6',1,2,1,.001,.01,'NoCytomix'],
              ['GENE7',1,2,1,.001,.01,'UNKNOWN']]
        factors={'TNF':{'gene_symbol':'TNF','full_panel_lr_eligible':'1'}}
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'test.xlsx'
            workbook(path,{'Stem':duplicate,'Goblet':duplicate,'TA':rows})
            audits,records,allowed=a.parse_workbook(path,factors,{'GENE1','GENE2','TNF','TNFRSF1A'})
        self.assertEqual(allowed,{'TA'})
        self.assertEqual(audits[2]['control_rows'],1)
        self.assertEqual(audits[2]['unknown_factor_rows'],1)
        self.assertEqual(audits[2]['invalid_or_missing_padj_rows'],1)
        result,bycell=a.factor_tables(factors,records,{'GENE1','GENE2','TNF','TNFRSF1A'},allowed)
        self.assertEqual(result[0]['n_strong_genes_in_panel'],1)
        self.assertEqual(result[0]['all_sheets_n_strong_in_panel'],2)
        self.assertEqual(next(x for x in bycell if x['receiver_sheet']=='TA')['n_down_in_panel'],1)
        self.assertEqual(result[0]['conservative_panel_candidate'],0)

    def test_receptor_requires_same_cell_all_subunits(self):
        import numpy as np
        positive={'A':np.array([1,0,1],bool),'B':np.array([0,1,1],bool)}
        self.assertEqual(a.complete_positive(positive,[('A','B')]).tolist(),[False,False,True])
        self.assertFalse(a.complete_positive(positive,[('A','MISSING')]).any())

    def test_malformed_workbook_fails(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'bad.xlsx'
            workbook(path,{'Bad':[['gene','pvalue'],['A',.01]]})
            with self.assertRaises(ValueError):a.parse_workbook(path,{'TNF':{}},set())

    def test_scientific_notation_text(self):
        self.assertTrue(a.significant(a.fnumber('4.52857414874191e-06'),a.fnumber('-1.01')))
        self.assertTrue(math.isnan(a.fnumber('NA')))

    def test_controls_and_real_gene_names(self):
        self.assertTrue(a.technical_target('Negative1'))
        self.assertTrue(a.technical_target('SystemControl324'))
        self.assertFalse(a.technical_target('NEGR1'))
        self.assertFalse(a.technical_target('CTRL'))

    def test_count_scan_qc_and_receptor_target_codetection(self):
        import gzip
        from unittest.mock import patch
        import specimens
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            (root/'counts').mkdir()
            header=['fov','cell_ID','cell','IFNG','IFNGR1','IFNGR2',
                    'G1','G2','G3','G4','G5','SystemControl1']
            data=[
                [1,1,'a',60,1,0,1,1,1,1,1,0],
                [1,2,'b',60,0,1,1,1,1,1,1,0],
                [1,3,'c',60,1,1,1,1,1,1,1,0],
                [1,4,'d',0,0,0,0,0,0,0,0,100],
                [2,5,'e',60,1,1,1,1,1,1,1,0]]
            with gzip.open(root/'counts/test_exprMat_file.csv.gz','wt',newline='') as f:
                w=csv.writer(f);w.writerow(header);w.writerows(data)
            rows=[dict(factor='IFNG',ligand_gene='IFNG',conservative_panel_candidate=1,
                       canonical_complexes_in_panel='IFNGR1_IFNGR2')]
            records=[dict(factor='IFNG',sheet='TA',strong=1,in_panel=1,gene=f'G{i}') for i in range(1,6)]
            with patch.object(specimens,'SPECIMENS',{'test':[('GEO','healthy',1,1)]}):
                gs,ss,cs,ps,qc=a.scan_expression(root,rows,records,{'TA'},50,2)
            self.assertEqual(qc,{'GEO':3})
            self.assertEqual(ss[0]['ligand_positive_cells'],3)
            self.assertEqual(ss[0]['canonical_receptor_positive_cells'],1)
            self.assertEqual(cs[0]['all_subunits_positive_cells'],1)
            self.assertEqual(ps[0]['cells_with_5_targets_detected'],3)
            self.assertEqual(ps[0]['canonical_receptor_and_5_targets_same_cell'],1)

    def test_empty_csv_has_schema(self):
        with tempfile.TemporaryDirectory() as d:
            m=a.save_table(Path(d),'empty.csv',[],['factor','gene'])
            self.assertEqual(m['rows'],0)
            self.assertGreater(m['bytes'],0)

if __name__=='__main__':unittest.main()
