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

        # PI-supplied provisional threshold in the current primary PCA.
        valid=int(row['group'])>=0
        hit=valid and float(row['pc1'])>3.1
        row['dust_pc1_rule_applicable']=valid
        row['dust_pc1_rule_match']=hit
        row['physical_label']='Dust impact' if hit or ref else 'Unclassified'
        row['physical_label_status']='user-identified reference' if ref else 'provisional' if hit else 'unclassified'
        row['physical_label_basis']='PI reference event' if ref else 'PI rule: PC1 > 3.1' if hit else ''
