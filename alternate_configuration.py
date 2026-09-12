"""Separate morphology and count/voltage exploration for V5 + SCM4 records."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy import signal,stats
from scipy.spatial.distance import cdist,pdist,squareform
from scipy.cluster.hierarchy import linkage,cut_tree
from cdflib import CDF
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_tds import P,CHANNELS,waveform,describe,dump
from build_atlas import simple_fit,savecsv,plot_event
from analyze_tds import robust_fit,silhouettes
from electron_voltage import norm,linear_lags

def run(args):
 rows=json.loads((args.cache/'final_records.json').read_text());prev=json.loads((args.cache/'previews.json').read_text());ids=[i for i,r in enumerate(rows) if r.get('analog_configuration')=='alternate_configuration'];fitids=[i for i in ids if not rows[i]['Burst_Saturation_Flag']];names=['V2','V5','V1V2','V3V4','SCM4'];XX=[];keys=None;amplitudes=[];corr=[];last=None
 for i in ids:
  r=rows[i]
  if last!=r['file']:c=CDF(args.data/r['file']);last=r['file']
  x=waveform(c,r['record'],r['N_Samples_per_Channel']);fe,am,aux,pv=describe(x,r['Sample_Speed'])
  keys=list(fe);XX.append(list(fe.values()));amplitudes.append([aux[ch+'__rms'] for ch in CHANNELS[:4]])
  if r['Burst_Saturation_Flag'] or r['counts_total']<32:continue
  centered=x[:4]-np.median(x[:4],axis=1,keepdims=True);env=np.sqrt(np.mean(centered.reshape(4,256,128)**2,axis=2));means=centered.reshape(4,256,128).mean(axis=2)
  raw=np.array(prev[i]['counts']['bins']);res=np.array(prev[i]['counts']['phase_residual']);lags=np.arange(-16,17);a=np.r_[env,means];rc=linear_lags(a,res,lags);hit=np.unravel_index(np.argmax(abs(rc)),rc.shape)
  item={k:r[k] for k in ['event','file','record','day','utc','counts_total','SWEAP_Start','solar_distance_Rs']};item.update(peak_channel=names[hit[0]%4],peak_measure='envelope' if hit[0]<4 else 'signed_mean',peak_residual_r=float(rc[hit]),counts_lag_ms=float(lags[hit[1]]*128/r['Sample_Speed']*1000))
  for ch,name in enumerate(names[:4]):
   item[name+'_raw_envelope_r0']=float(norm(env)[ch]@norm(raw));item[name+'_residual_envelope_r0']=float(rc[ch,16]);item[name+'_residual_signed_r0']=float(rc[ch+4,16])
  corr.append(item)
 X=np.array(XX);fit=np.array([ids.index(i) for i in fitids])
 try:
  lab,si,ag,Z,diag=simple_fit(X[fit],maxk=4);diag['accepted_partition']=True
 except ValueError as exc:
  if 'max() iterable argument is empty' not in str(exc):raise
  mid,scale=robust_fit(X[fit]);Z=np.clip((X[fit]-mid)/scale,-5,5)/np.sqrt(X.shape[1]);tree=linkage(Z,method='ward');D=squareform(pdist(Z));candidates=[]
  for k in range(2,5):
   ll=cut_tree(tree,n_clusters=k).ravel();candidates.append({'k':k,'sizes':np.bincount(ll).tolist(),'silhouette':float(silhouettes(D,ll).mean())})
  lab=np.zeros(len(fit),int);si=np.zeros(len(fit));ag=np.zeros(len(fit));diag={'k':1,'accepted_partition':False,'candidates':candidates,'silhouette':None,'resample_ARI_median':None,'interpretation':'No candidate partition met the minimum 20 records per group. A0 is an unpartitioned browsing pool, not a discovered class.'}
 for r in rows:r['alternate_group']=-1;r['alternate_status']='not applicable'
 groups=[]
 for j,i in enumerate(fitids):rows[i].update(alternate_group=int(lab[j]),alternate_status=('core' if si[j]>=.05 and ag[j]>=.8 else 'ambiguous') if diag['accepted_partition'] else 'unpartitioned')
 for g in range(diag['k']):
  loc=np.flatnonzero(lab==g);med=fitids[loc[np.argmin(cdist(Z[loc],Z[loc]).sum(axis=1))]];groups.append({'group':g,'n':len(loc),'core':sum(rows[fitids[j]]['alternate_status']=='core' for j in loc),'medoid':med})
 actualkeys=[k.replace('SCM5','SCM4').replace('V4__','V5__').replace('V2_V4','V2_V5').replace('V4_V1V2','V5_V1V2').replace('V4_V3V4','V5_V3V4').replace('V4_SCM4','V5_SCM4') for k in keys]
 # Replace names token-wise to avoid corrupting V3V4 feature prefixes.
 actualkeys=[]
 for key in keys:
  left,right=key.split('__');parts=left.split('_');left='_'.join({'V4':'V5','SCM5':'SCM4'}.get(v,v) for v in parts);actualkeys.append(left+'__'+right)
 amps=np.array(amplitudes)[fit];total=np.array([rows[i]['counts_total'] for i in fitids]);days=np.array([rows[i]['day'] for i in fitids]);across=[]
 for j,name in enumerate(names[:4]):
  across.append({'channel':name,'pooled_spearman':float(stats.spearmanr(amps[:,j],total).statistic),'by_day':{d:{'n':int(np.sum(days==d)),'spearman':float(stats.spearmanr(amps[days==d,j],total[days==d]).statistic) if np.sum(days==d)>=10 else None} for d in np.unique(days)},'median_residual_envelope_r0':float(np.median([r[name+'_residual_envelope_r0'] for r in corr]))})
 diag.update(channels=names,n_records=len(ids),n_fit=len(fitids),groups=groups,fit_indices=fitids,features=actualkeys,method='Separate Ward fit with robust scaling and equal per-feature weights, k=2..4 with >=20 records per group. Forty record/feature resamples only if a partition qualifies. A numbering is independent of G.',correlations={'n':len(corr),'across':across,'caveat':'Small, heterogeneous sample; different SWEAP status and distance from reference configuration. Pooled correlations do not isolate within-burst coupling.'})
 dump(args.out/'alternate_configuration_diagnostics.json',diag);savecsv(args.out/'alternate_configuration_features.csv',[dict(event=i,**dict(zip(actualkeys,X[j].tolist()))) for j,i in enumerate(ids)]);savecsv(args.out/'alternate_voltage_correlations.csv',corr);dump(args.cache/'final_records.json',rows)
 for g in groups:plot_event(args,rows[g['medoid']],args.figures/'examples'/f"alternate_event_{g['medoid']:04d}.png")
 fig,axs=plt.subplots(1,2,figsize=(12,5),layout='constrained');xx=np.arange(4)
 axs[0].bar(xx,[a['pooled_spearman'] for a in across]);axs[0].set_xticks(xx,names[:4]);axs[0].set(ylabel='Spearman correlation',title=f'Burst counts vs voltage RMS · n={len(fitids)}')
 for j,name in enumerate(names[:4]):axs[1].hist([r[name+'_residual_envelope_r0'] for r in corr],bins=np.linspace(-1,1,26),histtype='step',label=name)
 axs[1].set(xlabel='Phase-residual / voltage-envelope r',ylabel='Bursts',title=f'Within-burst, zero lag · n={len(corr)}');axs[1].legend();fig.suptitle('Alternate analog configuration · V5 / SCM4\nSeparate from V4 / SCM5 comparisons; scan baseline conditioned on SWEAP status')
 fig.savefig(args.figures/'alternate_voltage_correlations.png',dpi=150);plt.close(fig)
 lines=['# Alternate V5 / SCM4 configuration','',f'{len(ids)} records, {len(fitids)} unsaturated records fitted separately. Channels: '+', '.join(names)+'.', '',diag['method'],'',f"Selected {diag['k']} groups; silhouette {diag['silhouette']:.3f}; median resampling ARI {diag['resample_ARI_median']:.3f}." if diag['accepted_partition'] else diag['interpretation'],'','| A group | Records | Core | Medoid event |','|---|---:|---:|---:|']+[f"| A{g['group']} | {g['n']} | {g['core']} | {g['medoid']} |" for g in groups]
 lines+=['','## Electron–voltage correlations','',f"Within-burst analysis uses {len(corr)} unsaturated bursts with ≥32 counts. These coefficients use RMS envelopes in 128-sample bins and the sweep-phase baseline learned on other days with the same SWEAP status. No physical labels or calibration corrections are applied.",'','| Voltage | Pooled total / RMS Spearman | Median phase-residual envelope r at zero lag |','|---|---:|---:|']+[f"| {a['channel']} | {a['pooled_spearman']:.3f} | {a['median_residual_envelope_r0']:.3f} |" for a in across]
 lines+=['','Pooled associations are sensitive to date, distance and instrument state. Daily coefficients and per-event lag correlations are in the diagnostic JSON and alternate_voltage_correlations.csv. These groups are exploratory: a small sample cannot establish stable physical populations. The primary G fit excludes this configuration; the separate A fit retains it without treating V5 as V4 or SCM4 as SCM5.']
 (args.out/'alternate_configuration_report.md').write_text('\n'.join(lines));print(json.dumps(diag),flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser()
 for name in ['data','cache','out','figures']:ap.add_argument('--'+name,type=Path,required=True)
 run(ap.parse_args())
