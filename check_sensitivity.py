"""Persistent-line spectral sensitivity and cross-day SWEAP phase conditioning."""
import argparse, json, csv
from pathlib import Path
import numpy as np
from scipy import signal, ndimage
from scipy.spatial import cKDTree
from cdflib import CDF
from analyze_tds import P,CHANNELS,waveform,robust_fit,transform,labels_for,ari,dump
from build_atlas import simple_fit,savecsv

def run(args):
    rows=json.loads((args.cache/'final_records.json').read_text());prev=json.loads((args.cache/'previews.json').read_text())
    model=json.loads((args.out/'model_diagnostics.json').read_text());z=np.load(args.cache/'features.npz');X=z['X'];keys=z['keys'].tolist()
    ix=np.array(model['fit_indices']);labs=np.array([r['group'] for r in rows])[ix]
    # Full-resolution Welch spectra, one record at a time, cached for this diagnostic.
    if (args.cache/'spectra.npy').exists() and np.load(args.cache/'spectra.npy',mmap_mode='r').shape[0]==len(rows):spec=np.load(args.cache/'spectra.npy')
    else:
        spec=[];last=None
        if (args.cache/'spectra.npy').exists():spec=list(np.load(args.cache/'spectra.npy'))
        if len(spec)>len(rows):spec=[]
        for r in rows[len(spec):]:
            if last!=r['file']:c=CDF(args.data/r['file']);last=r['file']
            x=waveform(c,r['record'],r['N_Samples_per_Channel'])
            f,ps=signal.welch(x,r['Sample_Speed'],nperseg=4096,noverlap=2048,axis=1)
            spec.append(ps[:,1:].astype(np.float32))
        spec=np.array(spec);np.save(args.cache/'spectra.npy',spec)
    fs=rows[0]['Sample_Speed'];freq=np.arange(1,2049)*fs/4096
    logps=10*np.log10(np.maximum(spec,1e-30));background=ndimage.median_filter(logps,size=(1,1,31),mode='nearest')
    excess=logps-background;persistence=np.mean(excess[ix]>8,axis=0);medianex=np.median(excess[ix],axis=0)
    masks=[];lines=[]
    for ch in range(5):
        peaks,_=signal.find_peaks(medianex[ch],height=8,distance=3)
        peaks=[p for p in peaks if persistence[ch,p]>=.6 and p>=16]
        mask=np.zeros(len(freq),bool)
        for p in peaks:mask[max(0,p-1):p+2]=True
        masks.append(mask);lines.append({'channel':CHANNELS[ch],'frequencies_hz':freq[peaks].tolist(),
             'median_removed_power_fraction':float(np.median(spec[ix,ch][:,mask].sum(axis=1)/spec[ix,ch].sum(axis=1)))})
    newX=X.copy()
    for j in ix:
        for ch in range(5):
            keep=~masks[ch];f=freq[keep];p=spec[j,ch,keep].astype(float);p/=p.sum()
            fq=np.interp([.05,.5,.95],np.cumsum(p),f)
            vals={'spectral_entropy':-np.sum(p*np.log(np.maximum(p,1e-30)))/np.log(len(p)),
                  'log_f50':np.log10(fq[1]),'log_bandwidth':np.log10(max(fq[2]-fq[0],1)),
                  'spectral_peak3':np.sort(p)[-3:].sum()}
            for k,v in vals.items():newX[j,keys.index(CHANNELS[ch]+'__'+k)]=v
    mm,ss=robust_fit(newX[ix]);zz,_=transform(newX[ix],mm,ss,keys);ll=labels_for(zz,model['k'])
    line_result={'method':'Mask +/-1 Welch bin around peaks with local excess >8 dB in >=60% of unsaturated records; 31-bin local median, peaks above 7.5 kHz. Recompute four spectral descriptors per analog channel; retain all temporal and coupling descriptors.',
                 'channels':lines,'fixed_k_ARI_vs_raw':ari(labs,ll),
                 'caveat':'Sensitivity test only. Persistent lines are not proven instrumental, and masking can remove real signals. No time-domain filtering.'}
    for r in rows:r['line_mask_group']=-1
    for j,i in enumerate(ix):rows[i]['line_mask_group']=int(ll[j])
    # Phase profile is learned from binned counts on all other days, without physics labels.
    bins=np.array([v['counts']['bins'] for v in prev],float);phase=np.array([r['SWEAP_Start'] for r in rows])
    day=np.array([r['day'] for r in rows]);dt=128/fs
    samplephase=phase[:,None]+(np.arange(256)+.5)*dt
    edges=np.linspace(0,float(samplephase.max())+1e-8,601);center=(edges[:-1]+edges[1:])/2
    baseline=np.zeros_like(bins);models={}
    status=np.array([str(r['SWEAP_Status']) for r in rows])
    for d in np.unique(day):
      for state in np.unique(status[day==d]):
        target=(day==d)&(status==state);source=(day!=d)&(status==state)
        training_n=int(np.sum(source));fallback=training_n<20
        if fallback:source=day!=d
        xp=samplephase[source].ravel();yp=bins[source].ravel();bi=np.clip(np.digitize(xp,edges)-1,0,len(center)-1)
        profile=np.array([np.median(yp[bi==i]) if np.any(bi==i) else np.nan for i in range(len(center))])
        valid=np.isfinite(profile);profile=np.interp(center,center[valid],profile[valid])
        profile=ndimage.median_filter(profile,size=3,mode='nearest')
        baseline[target]=np.interp(samplephase[target],center,profile)
        for i in np.flatnonzero(target):rows[i]['counts_phase_baseline_status_matched']=not fallback;rows[i]['counts_phase_training_n']=training_n
        models[d+'_'+state]={'trained_on':sorted(set(day[source])),'SWEAP_Status':state,'cross_status_fallback':fallback,'matching_status_training_records':training_n,'minimum_training_records':20,'phase_centers_s':center.tolist(),'expected_counts_per_bin':profile.tolist()}
    predicted=baseline.sum(axis=1);truth=bins.sum(axis=1)
    baseline_valid=np.array([r['counts_phase_baseline_status_matched'] for r in rows]);yy=np.log1p(truth[baseline_valid]);yp=np.log1p(predicted[baseline_valid]);r2=float(1-np.sum((yy-yp)**2)/np.sum((yy-yy.mean())**2))
    usable=np.array([r['counts_usable_shape'] and r['counts_phase_baseline_status_matched'] for r in rows]);ci=np.flatnonzero(usable);resfeatures=[]
    for i in range(len(rows)):
        # Overall scale is retained separately; residual shape removes a per-burst gain.
        gain=(truth[i]+1)/(predicted[i]+1);expected=baseline[i]*gain
        residual=(bins[i]-expected)/np.sqrt(expected+1)
        centered=residual-residual.mean();sd=max(float(np.std(centered)),1e-12)
        rr=centered/sd;f,ps=signal.welch(rr,fs/128,nperseg=128,noverlap=64);p=ps[1:]/max(ps[1:].sum(),1e-30)
        features=[np.log10(max(np.sqrt(np.mean(residual**2)),1e-8)),np.log10(max(np.max(np.abs(rr)),1e-8)),
                  np.mean(rr**3)/(max(np.mean(rr**4),1e-12)**.75),np.corrcoef(rr[:-1],rr[1:])[0,1] if np.std(rr)>0 else 0,
                  -np.sum(p*np.log(np.maximum(p,1e-30)))/np.log(len(p)),np.sort(rr**2)[-26:].sum()/max(np.sum(rr**2),1e-20)]
        resfeatures.append(features);rows[i]['counts_phase_gain']=float(gain);rows[i]['counts_phase_expected_total']=float(predicted[i])
        rows[i]['counts_phase_residual_rms']=float(np.sqrt(np.mean(residual**2)))
        prev[i]['counts']['phase_expected']=np.round(expected,3).tolist();prev[i]['counts']['phase_residual']=np.round(residual,3).tolist()
    RF=np.array(resfeatures);cl,si,ag,zr,rd=simple_fit(RF[ci],maxk=6)
    for r in rows:r['residual_count_group']=-1;r['residual_count_status']='insufficient counts' if r['counts_phase_baseline_status_matched'] else 'insufficient matching-status training records'
    for j,i in enumerate(ci):rows[i].update(residual_count_group=int(cl[j]),residual_count_status='core' if si[j]>=.05 and ag[j]>=.8 else 'ambiguous')
    from scipy.spatial.distance import cdist
    rd['groups']=[]
    for g in range(rd['k']):
        sub=np.flatnonzero(cl==g);ids=ci[sub];med=int(ids[np.argmin(cdist(zr[sub],zr[sub]).sum(axis=1))])
        rd['groups'].append({'group':g,'n':len(ids),'medoid':med,'core':sum(rows[i]['residual_count_status']=='core' for i in ids),
                            'median_residual_rms':float(np.median([rows[i]['counts_phase_residual_rms'] for i in ids]))})
    rd.update(cross_day_log1p_total_R2=r2,baseline_valid_n=int(baseline_valid.sum()),baseline_cross_status_fallback_n=int((~baseline_valid).sum()),profiles=models,features=['log_residual_rms','log_crest','asymmetry','lag1_corr','spectral_entropy','top10pct_energy'],
        caveat='Cross-day empirical phase profile within SWEAP status, 600 bins; gain adjusted per event. Counts standardized by sqrt(expected+1), not a calibrated counting-noise model. Fine scan structure, timing drift, and rate changes can remain.')
    dump(args.out/'sensitivity_diagnostics.json',{'persistent_lines':line_result,'phase_conditioned_counts':rd})
    dump(args.cache/'final_records.json',rows);dump(args.cache/'previews.json',prev);savecsv(args.out/'event_catalog.csv',rows)
    savecsv(args.out/'phase_residual_features.csv',[dict(event=i,**dict(zip(rd['features'],RF[i].tolist()))) for i in range(len(rows))])
    print('LINE',json.dumps(line_result),flush=True)
    print('PHASE',json.dumps({k:v for k,v in rd.items() if k not in ['profiles']}),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);ap.add_argument('--cache',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);run(ap.parse_args())
