"""Combine the native-waveform morphology and harmonic-masked spectral features."""
from pathlib import Path
import json,csv
import numpy as np
from scipy.cluster.hierarchy import linkage,cut_tree
P=Path('outputs/tds_github_pages');O=Path('outputs/tds_waveform_atlas/combined_pca');O.mkdir(exist_ok=True)
rows=json.loads((P/'pooled_waveform_points.json').read_text());spectra={(r['event'],r['channel']):r for r in json.loads((P/'pooled_spectral_channels.json').read_text())}
wm=np.load('outputs/tds_waveform_atlas/pooled_waveform/model.npz');sm=np.load('outputs/tds_waveform_atlas/pooled_spectral/model.npz');keys=wm['features'].tolist();extra=json.loads(Path('outputs/tds_waveform_atlas/pooled_waveform/recovered_features.json').read_text())
caches={};indices={};source={}
for enc,cache in [('E23','work/tds'),('E24','work/tds_e24')]:
 caches[enc]=np.load(Path(cache)/'features.npz');indices[enc]=[caches[enc]['keys'].tolist().index(k) for k in keys]
 for i,r in enumerate(json.loads((Path(cache)/'records.json').read_text())):source[enc,r['file'],r['record']]=i
W=[];S=[];selected=[]
for r in rows:
 v=[spectra.get((r['event'],ch)) for ch in ['V2','V1V2','V3V4']]
 if not r['fit_included'] or not all(ch and ch['valid'] for ch in v):continue
 x=np.array(extra[str(r['event'])]) if str(r['event']) in extra else caches[r['encounter']]['X'][source[r['encounter'],r['file'],r['record']],indices[r['encounter']]]
 W.append(x);S.append([np.log10(max(ch[k],1e-30)) for ch in v for k in ['peak_hz','fractional_bandwidth','peak_power_fraction','peak_to_background']]);selected.append(r)
W=np.clip((np.array(W)-wm['median'])/wm['scale'],-5,5)*wm['weights'];S=np.clip((np.array(S)-sm['median'])/sm['scale'],-5,5)
assert np.isfinite(W).all() and np.isfinite(S).all()
# Preserve established per-feature scaling, then give each block half the
# total input variance on this common, complete-case population.
wo=W.mean(0);so=S.mean(0);W-=wo;S-=so;ws=np.sqrt(np.mean(np.sum(W*W,axis=1)));ss=np.sqrt(np.mean(np.sum(S*S,axis=1)));Z=np.c_[W/ws,S/ss]/np.sqrt(2)
u,s,v=np.linalg.svd(Z,full_matrices=False)
for component in v:
 if component[np.argmax(abs(component))]<0:component*=-1
scores=Z@v[:3].T;labels=cut_tree(linkage(Z,method='ward'),n_clusters=3).ravel()
assert np.max(abs((Z@v.T)@v-Z))<1e-10
block_variance=[float(np.var(Z[:,:27],axis=0).sum()),float(np.var(Z[:,27:],axis=0).sum())];assert np.allclose(block_variance,[.5,.5])
points=[]
for i,r in enumerate(selected):
 a=spectra[r['event'],'V1V2'];points.append(dict(event=r['event'],encounter=r['encounter'],utc=r['utc'],pc1=float(scores[i,0]),pc2=float(scores[i,1]),pc3=float(scores[i,2]),group=int(labels[i]),url='index.html#event='+str(r['event']),peak_hz=a['peak_hz'],df_f=a['fractional_bandwidth']))
ref=next(i for i,r in enumerate(selected) if r['event']==1786);d=np.linalg.norm(Z-Z[ref],axis=1);neighbors=[dict(event=selected[i]['event'],encounter=selected[i]['encounter'],distance=float(d[i])) for i in np.argsort(d) if i!=ref][:10]
summary=dict(total_events=len(rows),fitted_events=len(selected),fitted_by_encounter={enc:sum(r['encounter']==enc for r in selected) for enc in ['E23','E24']},explained_variance=(s[:3]**2/(s*s).sum()).tolist(),group_sizes=np.bincount(labels).tolist(),features=['waveform:'+k for k in keys]+['spectral:'+k for k in sm['features'].tolist()],input_block_variance=block_variance,pc_waveform_squared_loading_fraction=np.sum(v[:3,:27]**2,axis=1).tolist(),reference1786_neighbors=neighbors,method='39-feature combined PCA: 27 waveform descriptors plus 12 nonharmonic spectral descriptors, common V2/V1V2/V3V4 channels. Retain existing robust feature scaling/clipping and waveform family weights; center each block on the joint eligible population, normalize its total variance to 0.5, and concatenate. No missing-value imputation. Saturation and any invalid common-channel descriptor exclude an event from this fit. Equal waveform/spectral weighting is an explicit exploratory choice. Three Ward groups use all 39 combined features and are not physical classes. Original E23 dust thresholds are not applied to these axes.')
for name,data in [('points',points),('summary',summary)]:
 text=json.dumps(data,allow_nan=False,separators=(',',':'));(O/(name+'.json')).write_text(text);(P/('combined_pca_'+name+'.json')).write_text(text)
np.savez(O/'model.npz',waveform_median=wm['median'],waveform_scale=wm['scale'],waveform_weights=wm['weights'],spectral_median=sm['median'],spectral_scale=sm['scale'],waveform_origin=wo,spectral_origin=so,waveform_block_scale=ws,spectral_block_scale=ss,components=v,event_ids=np.array([r['event'] for r in selected]))
merged=json.loads((P/'merged_events.json').read_text());lookup={r['event']:r for r in points}
for r in merged:
 a=lookup.get(r['event']);r['combined_pca_fit_included']=a is not None
 for k in ['pc1','pc2','pc3','group']:r['combined_pca_'+k]=a[k] if a else None
(P/'merged_events.json').write_text(json.dumps(merged,allow_nan=False,separators=(',',':')))
with (P/'merged_events.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(merged[0]));w.writeheader();w.writerows(merged)
print(json.dumps(summary,indent=2))
