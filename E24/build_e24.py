"""Validate E24 extraction and project onto the frozen E23 basis."""
from pathlib import Path
import sys,json,csv,gzip,shutil
from collections import Counter
import numpy as np
from scipy.spatial.distance import cdist
sys.path.insert(0,'outputs/tds_waveform_atlas')
from analyze_tds import transform

O=Path('outputs/tds_waveform_atlas/E24');O.mkdir(exist_ok=True)
C=Path('work/tds_e24');OLD=Path('work/tds')
model=json.loads(Path('outputs/tds_waveform_atlas/model_diagnostics.json').read_text())
old=np.load(OLD/'features.npz');new=np.load(C/'features.npz');keys=old['keys'].tolist()
assert keys==new['keys'].tolist()
mid=np.array(model['scaling_median']);scale=np.array(model['scaling_iqr'])
z,_=transform(old['X'],mid,scale,keys);ix=np.array(model['fit_indices'])
origin=z[ix].mean(0);_,sv,v=np.linalg.svd(z[ix]-origin,full_matrices=False)
oldrows=json.loads((OLD/'final_records.json').read_text())
oldpc=np.array([[r['pc1'],r['pc2']] for r in oldrows]);error=float(np.max(abs((z-origin)@v[:2].T-oldpc)))
assert error<1e-8, error
rows=json.loads((C/'records.json').read_text());preview=json.loads((C/'previews.json').read_text())
zz,_=transform(new['X'],mid,scale,keys);pc=(zz-origin)@v[:3].T
centers=np.array(model['centroids']);dist=cdist(zz,centers)
# Stable ID allocation is persisted independently of arrival/file order.
registry=O/'event_ids.json'
ids=json.loads(registry.read_text()) if registry.exists() else {}
nextid=max([int(r['event']) for r in oldrows]+list(ids.values()))+1
for r in rows:
 key=r['file']+':'+str(r['record'])
 if key not in ids:ids[key]=nextid;nextid+=1
registry.write_text(json.dumps(ids,indent=2))
points=[]
for i,r in enumerate(rows):
 key=r['file']+':'+str(r['record']);r['local_event']=r['event'];r['event']=ids[key];r['encounter']='E24'
 valid=r['analog_configuration']=='reference' and np.isfinite(new['X'][i]).all()
 r['pca_projection_valid']=bool(valid)
 r['pc1']=float(pc[i,0]) if valid else None;r['pc2']=float(pc[i,1]) if valid else None
 r['reference_group']=int(dist[i].argmin()) if valid else None
 r['reference_centroid_distance']=float(dist[i].min()) if valid else None
 r['dust_rule_match']=bool(valid and pc[i,0]>3.1)
 r['physical_label']='Dust candidate (E23 threshold transferred)' if r['dust_rule_match'] else 'Unclassified'
 points.append({'event':r['event'],'utc':r['utc'],'file':r['file'],'record':r['record'],'pc1':r['pc1'],'pc2':r['pc2'],'group':r['reference_group'],'dust':r['dust_rule_match'],'saturated':bool(r['Burst_Saturation_Flag']),'config':r['analog_configuration'],'channels':r['analog_channels'],'counts_present':r['counts_present'],'counts':r['counts_total'] if r['counts_present'] else None,'fs':r['Sample_Speed'],'n':r['N_Samples_per_Channel'],'sweap_status':r['SWEAP_Status']})
inventory=json.loads((C/'inventory.json').read_text())
for a in inventory:
 f=Path('/Users/bale/Dropbox/projects/PSP/Work/TDS/data/E24')/a['file'];assert f.stat().st_size==a['bytes']
snapshot=json.loads(Path('work/e24_transfer_snapshot.json').read_text())
for a in snapshot:
 f=Path('/Users/bale/Dropbox/projects/PSP/Work/TDS/data/E24')/a['file'];assert f.stat().st_size==a['bytes'] and f.stat().st_mtime_ns==a['mtime_ns'],'Source changed during extraction'
summary={'encounter':'E24','source':'Dropbox','files':len(inventory),'events':len(rows),'id_min':min(r['event'] for r in rows),'id_max':max(r['event'] for r in rows),'utc_min':min(r['utc'] for r in rows),'utc_max':max(r['utc'] for r in rows),'configurations':dict(Counter(r['analog_configuration'] for r in rows)),'sample_rates':dict(Counter(r['Sample_Speed'] for r in rows)),'sample_lengths':dict(Counter(r['N_Samples_per_Channel'] for r in rows)),'saturation_flagged':sum(bool(r['Burst_Saturation_Flag']) for r in rows),'missing_counts':sum(not r['counts_present'] for r in rows),'dust_threshold_candidates':sum(r['dust_rule_match'] for r in rows),'projected_events':sum(r['pca_projection_valid'] for r in rows),'nearest_E23_group_counts':dict(Counter(r['reference_group'] for r in rows if r['pca_projection_valid'])),'E23_basis_reproduction_max_error':error,'method':'E24 projected using frozen E23 scaling, feature weights, PCA basis and nearest E23 cluster centroid. No new cluster fit. PC1>3.1 is an unvalidated threshold transfer to E24. Alternate and incomplete configurations have no PCA assignment. Raw spectral power-supply lines still affect the original E23 feature space.','input_files_stable_during_extraction':True}
(O/'summary.json').write_text(json.dumps(summary,indent=2));(O/'source_inventory.json').write_text(json.dumps(inventory,indent=2))
with (O/'event_catalog.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
manifest=[]
for start in range(0,len(preview),256):
 name=f'previews_{start//256:03d}.json.gz';manifest.append(name)
 with gzip.open(O/name,'wt',encoding='utf8',compresslevel=3) as f:json.dump(preview[start:start+256],f,allow_nan=False,separators=(',',':'))
(O/'preview_manifest.json').write_text(json.dumps(manifest))
(O/'events.json').write_text(json.dumps(points,allow_nan=False,separators=(',',':')))
html=Path('work/e24_template.html').read_text();(O/'index.html').write_text(html)
print(json.dumps(summary,indent=2),flush=True)
