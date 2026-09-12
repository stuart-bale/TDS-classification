"""Apply user labels by immutable source file and record, independently of clustering."""
import json
from pathlib import Path
def apply_expert_labels(rows):
    refs=json.loads((Path(__file__).parent/'expert_annotations.json').read_text())
    lookup={(v['file'],v['record']):v for v in refs}
    for row in rows:
        ref=lookup.get((row['file'],int(row['record'])))
        row['expert_label']=ref['label'] if ref else ''
        row['expert_label_source']=ref['source'] if ref else ''
