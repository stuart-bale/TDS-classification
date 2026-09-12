"""Publish a readable offline methods page from the current diagnostics."""
import json
from pathlib import Path

def build_methods(out):
 out=Path(out);inventory=json.loads((out/'source_inventory.json').read_text());model=json.loads((out/'model_diagnostics.json').read_text());phase=json.loads((out/'sensitivity_diagnostics.json').read_text())['phase_conditioned_counts']
 values={'__N__':f"{sum(f['records'] for f in inventory):,}",'__FILES__':str(len(inventory)),'__K__':str(model['k']),'__R2__':f"{phase['cross_day_log1p_total_R2']:.3f}",'__FALLBACK__':str(phase['baseline_cross_status_fallback_n'])}
 text=(Path(__file__).parent/'methods_template.html').read_text()
 for key,value in values.items():text=text.replace(key,value)
 (out/'methods.html').write_text(text)

if __name__=='__main__':build_methods(Path(__file__).parent)
