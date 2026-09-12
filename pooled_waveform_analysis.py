"""Joint encounter waveform morphology PCA; cached native-sample descriptors."""
from pathlib import Path
import json,sys
import numpy as np
from scipy.cluster.hierarchy import linkage,cut_tree
from cdflib import CDF
sys.path.insert(0,'outputs/tds_waveform_atlas')
from analyze_tds import describe,waveform,robust_fit
P=Path('outputs/tds_github_pages'); O=Path('outputs/tds_waveform_atlas/pooled_waveform');O.mkdir(exist_ok=True)
chs=['V2','V1V2','V3V4']
metrics=['log_crest','log_kurtosis','signed_asymmetry','energy_width','energy_asymmetry','envelope_occupancy','top1pct_energy']
keys=[ch+'__'+m for ch in chs for m in metrics]+[f'pair_{a}_{b}__{m}' for i,a in enumerate(chs) for b in chs[i+1:] for m in ['abs_corr','env_corr']]
rows=[];matrix=[];last=None;c=None
extra_path=O/'recovered_features.json';extra=json.loads(extra_path.read_text()) if extra_path.exists() else {}
ids=json.loads((P/'E24/event_ids.json').read_text())
oldmerged={r['event']:r for r in json.loads((P/'merged_events.json').read_text())}
e23meta={r['event']:r for r in json.loads((P/'atlas_metadata.json').read_text())['rows']}
for enc,cache in [('E23','work/tds'),('E24','work/tds_e24')]:
 raw=json.loads((Path(cache)/('final_records.json' if enc=='E23' else 'records.json')).read_text());z=np.load(Path(cache)/'features.npz');indices=[z['keys'].tolist().index(k) for k in keys]
 for i,r in enumerate(raw):
  if enc=='E23':
   if i%50==0:batch=json.loads((P/f'preview_chunks/{i//50}.json').read_text())
   r['analog_channels']=batch[i%50]['channels'];r['counts_present']=e23meta[r['event']]['counts_present']
   r['analog_configuration']='reference' if r['analog_channels']==['V2','V4','V1V2','V3V4','SCM5'] else 'alternate_configuration'
  event=r['event'] if enc=='E23' else ids[r['file']+':'+str(r['record'])];v=z['X'][i,indices];present=all(ch in r['analog_channels'] for ch in chs)
  if present and not np.isfinite(v).all():
   if str(event) in extra:v=np.array(extra[str(event)])
   else:
    file=(enc,r['file'])
    if file!=last:c=CDF(Path('/Users/bale/Dropbox/projects/PSP/Work/TDS/data')/enc/r['file']);last=file
    fe,_,_,_=describe(waveform(c,r['record'],r['N_Samples_per_Channel']),r['Sample_Speed']);v=np.array([fe[k] for k in keys]);extra[str(event)]=v.tolist()
  finite=bool(present and np.isfinite(v).all());saturated=bool(r['Burst_Saturation_Flag']);fit=finite and not saturated
  matrix.append(v);rows.append(dict(event=event,encounter=enc,utc=r['utc'],file=r['file'],record=r['record'],channels=r['analog_channels'],config=r['analog_configuration'],n=r['N_Samples_per_Channel'],fs=r['Sample_Speed'],counts_present=r['counts_present'],counts=r['counts_total'] if r['counts_present'] else None,sweap_status=r['SWEAP_Status'],saturated=saturated,dust=oldmerged[event]['dust'],preview_url=('' if enc=='E23' else 'E24/')+f'preview_chunks/{i//50}.json',preview_offset=i%50,fit_included=fit,projection_valid=finite,exclusion_reason=None if fit else 'saturation' if finite else 'missing or invalid common-channel features'))
  if i%1000==0:print(enc,i,flush=True)
extra_path.write_text(json.dumps(extra))
X=np.array(matrix);fit=np.array([r['fit_included'] for r in rows]);valid=np.array([r['projection_valid'] for r in rows]);mid,scale=robust_fit(X[fit]);weights=np.array([1/np.sqrt(6) if k.startswith('pair_') else 1/np.sqrt(21) for k in keys]);Z=np.clip((X-mid)/scale,-5,5)*weights;origin=Z[fit].mean(0);_,s,v=np.linalg.svd(Z[fit]-origin,full_matrices=False)
# Deterministic sign convention: largest absolute loading is positive.
for component in v:
 if component[np.argmax(abs(component))]<0:component*=-1
pc=(Z[valid]-origin)@v[:3].T
for r,coords in zip([r for r in rows if r['projection_valid']],pc):r.update(dict(zip(['pc1','pc2','pc3'],map(float,coords))))
for r in rows:
 for k in ['pc1','pc2','pc3']:r.setdefault(k,None)
labels=cut_tree(linkage(Z[fit],method='ward'),n_clusters=3).ravel();centers=np.array([Z[fit][labels==k].mean(0) for k in range(3)])
for r,label in zip([r for r in rows if r['fit_included']],labels):r['group']=int(label)
for r in rows:r.setdefault('group',None)
summary=dict(total_events=len(rows),fitted_events=int(fit.sum()),fitted_by_encounter={enc:sum(r['fit_included'] and r['encounter']==enc for r in rows) for enc in ['E23','E24']},projected_events=int(valid.sum()),features=keys,explained_variance=(s[:3]**2/(s*s).sum()).tolist(),group_sizes=np.bincount(labels).tolist(),method='Joint E23/E24 fit of 21 temporal shape descriptors and 6 common-channel correlations. V2, V1V2, V3V4 only; median removal and per-channel RMS normalization. Median/IQR scaling and clipping at 5 IQR, equal temporal/coupling family weight. No frequency descriptors, amplitudes, digital counts, SCM, or configuration flags enter this fit. Saturated events are projected for inspection but excluded from fitting and groups. Three Ward groups are exploratory, not physical labels. Original E23 PC1 dust annotations are retained; the 3.1 threshold is not applied to new axes.')
np.savez(O/'model.npz',median=mid,scale=scale,weights=weights,origin=origin,components=v,features=np.array(keys))
for name,data in [('points',rows),('summary',summary)]:
 text=json.dumps(data,allow_nan=False,separators=(',',':'));(O/(name+'.json')).write_text(text);(P/('pooled_waveform_'+name+'.json')).write_text(text)
for r in rows:
 m=oldmerged[r['event']];m['pooled_waveform_fit_included']=r['fit_included']
 for key in ['pc1','pc2','pc3','group']:m['pooled_waveform_'+key]=r[key]
(P/'merged_events.json').write_text(json.dumps(list(oldmerged.values()),allow_nan=False,separators=(',',':')))
import csv
with (P/'merged_events.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(next(iter(oldmerged.values()))));w.writeheader();w.writerows(oldmerged.values())
print(json.dumps(summary,indent=2),flush=True)
