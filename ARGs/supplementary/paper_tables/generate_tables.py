import csv,re
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from xml.sax.saxutils import escape
OUT=Path(__file__).resolve().parent
source=OUT.parents[1]/'tables/clinical_ARG_taxonomic_context_summary.csv'
raw=list(csv.DictReader(source.open()))
fields=['patient_id','follow_up_week','ARG','tools_detecting_ARG','n_tools','n_ARG_contigs','ARG_contigs','taxonomic_context','taxonomic_rank','lineage','min_coverage_percentage','max_sequence_identity']
# Keep manuscript visit labels while analysis tables store numeric visits.
s4=[[f"W{int(r[k])}" if k == "follow_up_week" else r[k] for k in fields] for r in raw]
assert len(s4)==19
CURATION = OUT / 'curation'


def read_curated_csv(filename, required_columns):
    path = CURATION / filename
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        missing = set(required_columns) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f'{path}: missing columns {sorted(missing)}')
        rows = list(reader)
    if not rows:
        raise ValueError(f'{path}: no curated records')
    return rows


contig_rows = read_curated_csv('curated_contigs.csv', [
    'sample', 'contig_id', 'contig_length_bp', 'ARG', 'ARG_start',
    'ARG_end', 'ARG_strand', 'MMseqs2_assignment', 'MMseqs2_rank'
])
cases = [
    (r['sample'], r['contig_id'], int(r['contig_length_bp']), r['ARG'],
     int(r['ARG_start']), int(r['ARG_end']), r['ARG_strand'],
     r['MMseqs2_assignment'], r['MMseqs2_rank'])
    for r in contig_rows
]
case_keys = {(r[0], r[1]) for r in cases}
if len({(r[0], r[1], r[3]) for r in cases}) != len(cases):
    raise ValueError('Duplicate sample-contig-ARG record in curated_contigs.csv')
contig_lengths = {}
for sample, cid, length, gene, start, end, strand, taxon, rank in cases:
    if not (1 <= start <= end <= length) or strand not in ('+', '-'):
        raise ValueError(f'Invalid ARG coordinates or strand: {sample}/{cid}')
    key = (sample, cid)
    if key in contig_lengths and contig_lengths[key] != length:
        raise ValueError(f'Conflicting contig lengths: {sample}/{cid}')
    contig_lengths[key] = length
query_rows = read_curated_csv('blast_queries.csv', [
    'sample', 'contig_id', 'BLAST_query_region', 'query_start_in_contig',
    'query_end_in_contig', 'evidence_available', 'source_reports',
    'source_report_order'
])
queries = {}
report_records = []
query_keys = set()
for r in query_rows:
    key = (r['sample'], r['contig_id'])
    if key not in case_keys:
        raise ValueError(f'Query has no curated contig: {key}')
    start, end = int(r['query_start_in_contig']), int(r['query_end_in_contig'])
    if not (1 <= start <= end <= contig_lengths[key]):
        raise ValueError(f'Invalid query coordinates: {key}')
    query_key = (*key, r['BLAST_query_region'], start, end)
    if query_key in query_keys:
        raise ValueError(f'Duplicate BLAST query: {query_key}')
    query_keys.add(query_key)
    queries.setdefault(key, []).append((
        r['BLAST_query_region'], start, end, r['evidence_available']
    ))
    filenames = [x.strip() for x in r['source_reports'].split(';') if x.strip()]
    orders = [int(x.strip()) for x in r['source_report_order'].split(';') if x.strip()]
    if len(filenames) != len(orders):
        raise ValueError(f'Mismatched report names and order: {key}')
    for order, filename in zip(orders, filenames):
        if Path(filename).name != filename:
            raise ValueError(f'Report must be a filename within alignment_reports: {filename}')
        report_records.append((order, *key, filename))
if set(queries) != case_keys:
    raise ValueError('Every curated contig must have at least one BLAST query')
if len({r[0] for r in report_records}) != len(report_records):
    raise ValueError('Duplicate source_report_order in blast_queries.csv')
if len({r[3] for r in report_records}) != len(report_records):
    raise ValueError('Duplicate alignment report in blast_queries.csv')
ha=['sample','contig_id','contig_length_bp','ARG','ARG_start','ARG_end','ARG_strand','MMseqs2_assignment','MMseqs2_rank','BLAST_query_region','query_start_in_contig','query_end_in_contig','query_length_bp','evidence_available']
a=[]
for sample,cid,n,g,start,end,strand,tax,rank in cases:
 for region,qs,qe,status in queries[(sample,cid)]:
  a.append([sample,cid,n,g,start,end,strand,tax,rank,region,qs,qe,qe-qs+1,status])
# Detailed BLAST HSPs; coordinates are retained in the reported direction.
reports = [(sample, cid, filename)
           for _, sample, cid, filename in sorted(report_records)]
hb=['sample','contig_id','query_region','query_length_bp','reference_accession','reference_description','reference_length_bp','HSP_number','query_start','query_end','reference_start','reference_end','strand','alignment_length_bp','identical_bases','identity_percent','gap_positions','query_coverage_this_HSP_percent','bit_score','E_value','source_report']
b=[]
archive=OUT/'alignment_reports'
missing_reports = [filename for _, _, filename in reports if not (archive / filename).is_file()]
if missing_reports:
 raise FileNotFoundError(f'Missing BLAST reports in {archive}: {missing_reports}')
for sample,cid,filename in reports:
 p=archive/filename
 s=p.read_text()
 acc=re.search(r'Sequence ID:\s*(\S+) Length:\s*(\d+)',s)
 desc=re.search(r'^>(.*)',s,re.M).group(1)
 qlen=int(re.search(r'Query ID:.*?Length:\s*(\d+)',s).group(1))
 for number,ch in enumerate(re.split(r'(?=Range \d+:)',s)[1:],1):
  q=re.findall(r'^Query\s+(\d+)\s+([A-Za-z-]+)\s+(\d+)',ch,re.M);t=re.findall(r'^Sbjct\s+(\d+)\s+([A-Za-z-]+)\s+(\d+)',ch,re.M)
  if not q:continue
  ident=re.search(r'Identities:\s*(\d+)/(\d+)',ch);gaps=re.search(r'Gaps:\s*(\d+)/(\d+)',ch);strand=re.search(r'Strand:\s*([^\n]+)',ch).group(1).strip()
  score=re.search(r'Score:\s*([\d.]+)',ch).group(1);ev=re.search(r'Expect:\s*([^,\s]+)',ch).group(1)
  qs,qe=int(q[0][0]),int(q[-1][2]);ss,se=int(t[0][0]),int(t[-1][2]);num,den=map(int,ident.groups())
  b.append([sample,cid,'complete contig',qlen,acc.group(1),desc,int(acc.group(2)),number,qs,qe,ss,se,strand,den,num,round(100*num/den,5),int(gaps.group(1)),round(100*(abs(qe-qs)+1)/qlen,5),float(score),ev,filename])

# Reviewed biological findings added to the supplementary tables.
interpretation_rows = read_curated_csv('contig_interpretations.csv', [
    'sample', 'contig_id', 'BLAST_and_sequence_comparison_findings',
    'interpretation_and_limitations'
])
findings = {}
for r in interpretation_rows:
    key = (r['sample'], r['contig_id'])
    if key in findings:
        raise ValueError(f'Duplicate contig interpretation: {key}')
    findings[key] = (
        r['BLAST_and_sequence_comparison_findings'],
        r['interpretation_and_limitations']
    )
if set(findings) != case_keys:
    raise ValueError('Interpretations must match the curated sample-contig records')
ha.extend(['BLAST_and_sequence_comparison_findings','interpretation_and_limitations'])
for row in a:
 row.extend(findings[(row[0], row[1])])
hb.insert(6,'reference_context')
for row in b:
 description=row[5].lower()
 context='plasmid' if 'plasmid' in description else ('chromosome' if 'chromosome' in description else ('transposon / class 1 integron' if 'integron' in description or 'transposon' in description else 'not specified'))
 row.insert(6,context)

notes=[
['S4','Original MMseqs2 labels and lineages are retained; they are contig classifications, not confirmed hosts.'],
['S4','19 patient-week-gene combinations correspond to 18 sample-specific contigs. PAC3-W0 k141_1044 contains two panel genes.'],
['S5A','One row per ARG and BLAST query. PAC3-W0 appears twice because the same contig contains two ARGs.'],
['S5A','Left/right refer to original contig coordinates, not transcriptional upstream/downstream. Coordinates are 1-based and inclusive.'],
['S5A','Ten distinct contigs were selected for exploratory review; this is not a uniform analysis of all 18 contigs.'],
['S5B','One row per reported HSP, including secondary/overlapping HSPs. Do not sum HSP lengths to infer nonredundant coverage.'],
['S5B','Identity is calculated from exact identical-base/alignment-length counts. Query coverage is for this HSP only.'],
['S5B','Detailed alignment files were available for six contigs; original flanking-query Descriptions summaries are not substituted for HSP-level data.'],
['Methods','BLASTn; database core_nt. Record exact run dates, BLAST version/task and parameter settings from original searches before submission. These are not inferred here.'],
['Interpretation','A reference plasmid/chromosome label does not establish the location or host in the metagenomic sample.'],
['Interpretation','Discontinuous reference matches are not, by themselves, evidence of an assembly error or biological rearrangement.'],
['Selection','Reference alignments were selected during manual exploratory review; not a claim that all matching references were exhaustively evaluated.']]
sheets=[('S4_taxonomic_context',fields,s4),('S5A_contigs_queries',ha,a),('S5B_reference_alignments',hb,b),('Notes',['scope','note'],notes)]
for name,h,rows in sheets:
 with (OUT/(name+'.csv')).open('w',newline='',encoding='utf-8-sig') as f:
  w=csv.writer(f);w.writerow(h);w.writerows(rows)
# SpreadsheetML workbook, written without external dependencies.
def col(n):
 s=''
 while n:n,r=divmod(n-1,26);s=chr(65+r)+s
 return s
ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
files={}
files['[Content_Types].xml']='<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'+''.join(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' for i in range(1,5))+'</Types>'
files['_rels/.rels']='<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>'
files['xl/workbook.xml']=f'<workbook xmlns="{ns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'+''.join(f'<sheet name="{name}" sheetId="{i}" r:id="rId{i}"/>' for i,(name,h,r) in enumerate(sheets,1))+'</sheets></workbook>'
files['xl/_rels/workbook.xml.rels']='<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'+''.join(f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>' for i in range(1,5))+'<Relationship Id="rId5" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>'
files['xl/styles.xml']=f'<styleSheet xmlns="{ns}"><fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><color rgb="FFFFFFFF"/><sz val="11"/><name val="Calibri"/></font></fonts><fills count="3"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill><fill><patternFill patternType="solid"><fgColor rgb="FF173047"/><bgColor indexed="64"/></patternFill></fill></fills><borders count="1"><border/></borders><cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs><cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"><alignment vertical="top" wrapText="1"/></xf><xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0"><alignment vertical="top" wrapText="1"/></xf></cellXfs><cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles></styleSheet>'
for i,(name,h,rows) in enumerate(sheets,1):
 xml=[]
 for rn,row in enumerate([h]+rows,1):
  cells=[]
  for cn,v in enumerate(row,1):
   ref=f'{col(cn)}{rn}';style=1 if rn==1 else 0
   if isinstance(v,(int,float)):cells.append(f'<c r="{ref}" s="{style}"><v>{v}</v></c>')
   else:cells.append(f'<c r="{ref}" s="{style}" t="inlineStr"><is><t xml:space="preserve">{escape(str(v))}</t></is></c>')
  xml.append(f'<row r="{rn}" ht="{45 if rn==1 else 60}" customHeight="1">'+''.join(cells)+'</row>')
 cols=''.join(f'<col min="{j}" max="{j}" width="{55 if any(x in key for x in ["description","lineage","note","evidence","findings","limitations"]) else 23}" customWidth="1"/>' for j,key in enumerate(h,1))
 files[f'xl/worksheets/sheet{i}.xml']=f'<worksheet xmlns="{ns}"><sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews><cols>{cols}</cols><sheetData>'+''.join(xml)+f'</sheetData><autoFilter ref="A1:{col(len(h))}{len(rows)+1}"/></worksheet>'
import xml.etree.ElementTree as ET
for v in files.values():ET.fromstring(v)
with ZipFile(OUT/'Supplementary_Tables_S4_S5.xlsx','w',ZIP_DEFLATED) as z:
 for name,data in files.items():z.writestr(name,data)
assert len({(r[0],r[1]) for r in a})==10
print('S4:',len(s4),'rows; S5A:',len(a),'rows / 10 contigs; S5B:',len(b),'HSPs; workbook XML verified.')
