"""Unsupervised TDS morphology exploration. Reads CDFs; never modifies them.

Dependencies: numpy, scipy, matplotlib, cdflib. Python >=3.10.
Run: python analyze_tds.py --data /path/to/cdfs --cache /path/to/cache --out /path/to/output
"""
import argparse, csv, hashlib, json, sys, time
from pathlib import Path
import numpy as np
from scipy import signal, stats
from scipy.cluster.hierarchy import linkage, cut_tree
from scipy.spatial.distance import pdist, squareform, cdist
from scipy.optimize import linear_sum_assignment
from cdflib import CDF, cdfepoch

P = 'PSP_FLD_L2_TDS_WF_'
CHANNELS = ['V2', 'V4', 'V1V2', 'V3V4', 'SCM5']
SUFFIXES = [x+'_Engineering_'+('nT' if x == 'SCM5' else 'mV') for x in CHANNELS]
SEED = 20250322

def dump(path, obj):
    Path(path).write_text(json.dumps(obj, allow_nan=False, separators=(',', ':')))

def rounded(a, n=5):
    return np.round(np.asarray(a), n).tolist()

def channel_names(c, rec):
    result=[]
    for primary,alternate in [('V2',None),('V4','V5'),('V1V2',None),('V3V4',None),('SCM5','SCM4')]:
        chosen=None
        for ch in [primary,alternate]:
            if ch and int(c.varget(P+'Burst_Time_Series_'+ch+'_Exists_Flag',startrec=rec,endrec=rec).ravel()[0])==1:
                chosen=ch;break
        result.append(chosen)
    return result

def waveform(c, rec, n):
    # Missing slots are zeros solely for constructing previews. Such records are
    # excluded from analog/joint fits and all their analog features are NaN.
    vals=[]
    for ch in channel_names(c,rec):
        if ch is None:
            vals.append(np.zeros(n));continue
        var=P+'Burst_Time_Series_'+ch+'_Engineering_'+('nT' if ch.startswith('SCM') else 'mV')
        x=c.varget(var,startrec=rec,endrec=rec).ravel()[:n].astype(float)
        fill=c.varattsget(var).get('FILLVAL',-1e31)
        if len(x)!=n or np.any(~np.isfinite(x)) or np.any(x==fill):
            raise ValueError(f'Invalid waveform: record {rec}, {var}')
        vals.append(x)
    return np.array(vals)

def describe(x, fs):
    # Median removal only: preserve slow asymmetric morphology.
    x=x-np.median(x, axis=1, keepdims=True)
    rms=np.sqrt(np.mean(x*x, axis=1))
    z=x/np.maximum(rms[:,None], 1e-20)
    f, ps=signal.welch(x, fs, window='hann', nperseg=4096, noverlap=2048,
                       detrend='constant', scaling='density', axis=1)
    ps=ps[:,1:]; f=f[1:]
    prob=ps/np.maximum(ps.sum(axis=1,keepdims=True), 1e-30)
    feat={}; amps={}; auxiliary={}
    env=np.sqrt(np.mean(z.reshape(5,128,-1)**2,axis=2))
    t=np.arange(x.shape[1])/fs
    for j,ch in enumerate(CHANNELS):
        zz=z[j]; energy=zz*zz; total=energy.sum()
        q=np.interp([.05,.5,.95], np.cumsum(energy)/max(total,1e-30), t)
        p=prob[j]; c=np.cumsum(p)
        fq=np.interp([.05,.5,.95], c, f)
        ac=signal.fftconvolve(zz, zz[::-1], mode='full')[len(zz)-1:]
        ac=ac/max(ac[0],1e-30)
        crossings=np.flatnonzero(ac<=0)
        start=max(3,int(crossings[0]) if len(crossings) else 3)
        acmax=float(np.max(ac[start:4096])) if start<4096 else 0.
        features={
            'log_crest':np.log10(max(np.max(np.abs(zz)),1e-12)),
            'log_kurtosis':np.log10(max(np.mean(zz**4),1e-12)),
            'signed_asymmetry':np.mean(zz**3)/(max(np.mean(zz**4),1e-12)**.75),
            'energy_width':(q[2]-q[0])/(len(zz)/fs),
            'energy_asymmetry':(q[2]+q[0]-2*q[1])/max(q[2]-q[0],1/fs),
            'envelope_occupancy':np.mean(env[j]>.5*env[j].max()),
            'top1pct_energy':np.sort(energy)[-max(1,len(zz)//100):].sum()/max(total,1e-30),
            'zero_crossing':np.mean(zz[:-1]*zz[1:]<0),
            'autocorr_peak':acmax,
            'spectral_entropy':-np.sum(p*np.log(np.maximum(p,1e-30)))/np.log(len(p)),
            'log_f50':np.log10(max(fq[1],1)),
            'log_bandwidth':np.log10(max(fq[2]-fq[0],1)),
            'spectral_peak3':np.sort(p)[-3:].sum(),
        }
        feat.update({ch+'__'+k:float(v) for k,v in features.items()})
        amps[ch+'__log_rms']=float(np.log10(max(rms[j],1e-20)))
        auxiliary[ch+'__rms']=float(rms[j])
        auxiliary[ch+'__peak_hz']=float(f[np.argmax(p)])
        auxiliary[ch+'__f50_hz']=float(fq[1])
    cc=np.corrcoef(z); ec=np.corrcoef(env)
    for a in range(5):
        for b in range(a+1,5):
            feat[f'pair_{CHANNELS[a]}_{CHANNELS[b]}__abs_corr']=float(abs(cc[a,b]))
            feat[f'pair_{CHANNELS[a]}_{CHANNELS[b]}__env_corr']=float(ec[a,b])
    # Average over 15 Welch segments, not a one-segment coherence estimate.
    for a,b in [(0,1),(2,3),(2,4),(3,4)]:
        fc,co=signal.coherence(x[a],x[b],fs,window='hann',nperseg=4096,noverlap=2048)
        weight=np.sqrt(prob[a]*prob[b]); weight/=max(weight.sum(),1e-30)
        feat[f'pair_{CHANNELS[a]}_{CHANNELS[b]}__coherence']=float(np.sum(co[1:]*weight))
    # Diagnose fixed instrumental lines, not used as a clustering feature.
    line=np.zeros(len(f),dtype=bool)
    for hz in np.arange(150000,fs/2,50000):
        line |= abs(f-hz)<=fs/4096
    auxiliary['fixed_line_fraction']=float(np.mean(prob[:,line].sum(axis=1)))
    # Viewer shows envelopes over the entire record, plus real samples in a zoom.
    blocks=x.reshape(5,512,-1)
    center=int(np.argmax(np.mean(z[:4]**2,axis=0)))
    lo=int(np.clip(center-512,0,len(t)-1024))
    edges=np.unique(np.rint(np.geomspace(1,len(f),180)).astype(int))-1
    psdb=[]; pf=[]
    for left,right in zip(edges[:-1],edges[1:]):
        if right>left:
            pf.append(float(np.sqrt(f[left]*f[right-1])))
            psdb.append(10*np.log10(np.maximum(np.mean(ps[:,left:right],axis=1),1e-30)))
    preview={'lo':rounded(blocks.min(axis=2)), 'hi':rounded(blocks.max(axis=2)),
             'zoom':rounded(x[:,lo:lo+1024]), 'zoom_start':lo,
             'freq':rounded(pf,2), 'psd_db':rounded(np.array(psdb).T,3),
             'rms':rounded(rms,6)}
    return feat,amps,auxiliary,preview

def describe_counts(c, analog, fs):
    n=len(c); total=float(c.sum()); mu=float(c.mean()); variance=float(c.var(ddof=1))
    if np.any(c<0) or np.any(~np.isfinite(c)) or not np.all(c==np.rint(c)):
        raise ValueError('Invalid digital counts')
    bins=c.reshape(256,-1).sum(axis=1)
    features={}
    for b in [16,128,1024]:
        sums=c.reshape(-1,b).sum(axis=1)
        nullvar=b*variance*(n-b)/(n-1)
        features[f'counts__log_dispersion_{b}']=float(np.log10(max(sums.var(ddof=1)/max(nullvar,1e-20),1e-8))) if total else 0.
    p=bins/max(total,1);q=np.interp([.05,.5,.95],np.cumsum(p),np.arange(256)/256) if total else [0,0,0]
    features['counts__width90']=float(q[2]-q[0])
    features['counts__asymmetry']=float((q[2]+q[0]-2*q[1])/max(q[2]-q[0],1/256))
    features['counts__top10pct_fraction']=float(np.sort(bins)[-26:].sum()/max(total,1))
    features['counts__bin_entropy']=float(-np.sum(p*np.log(np.maximum(p,1e-30)))/np.log(256))
    f,ps=signal.welch(c,fs,nperseg=4096,noverlap=2048)
    pp=ps[1:]/max(ps[1:].sum(),1e-30)
    features['counts__spectral_entropy']=float(-np.sum(pp*np.log(np.maximum(pp,1e-30)))/np.log(len(pp)))
    features['counts__spectral_peak3']=float(np.sort(pp)[-3:].sum())
    features['counts__fano_sample']=float(variance/max(mu,1e-20))
    # Cross-channel association uses coarsened count totals and analog RMS envelopes.
    analog=analog-np.median(analog,axis=1,keepdims=True)
    env=np.sqrt(np.mean(analog.reshape(5,256,-1)**2,axis=2))
    normenv=env-env.mean(axis=1,keepdims=True);normenv/=np.maximum(np.linalg.norm(normenv,axis=1,keepdims=True),1e-20)
    bc=bins-bins.mean();bc/=max(np.linalg.norm(bc),1e-20)
    corr=np.fft.irfft(np.fft.rfft(normenv,axis=1)*np.conj(np.fft.rfft(bc)),n=256,axis=1)
    for j,ch in enumerate(CHANNELS):features[f'count_pair_{ch}__envelope_corr']=float(corr[j,0])
    lagidx=np.arange(-16,17)%256;obs=float(np.max(np.abs(corr[:,lagidx])))
    # Circular-shift reference preserves the count envelope's autocorrelation.
    # Exploratory only: a scanning/nonstationary instrument can violate this null.
    null=np.array([np.max(np.abs(corr[:,(lagidx+s)%256])) for s in range(33,224)])
    pshift=float((1+np.sum(null>=obs))/(1+len(null))) if total else 1.
    aux={'counts_total':int(total),'counts_rate_cps':total/(n/fs),
         'counts_nonzero_fraction':float(np.mean(c>0)), 'counts_max_sample':int(c.max()),
         'counts_envelope_max_abs_corr':obs,'counts_shift_p_exploratory':pshift,
         'counts_usable_shape':bool(total>=32)}
    preview={'bins':bins.astype(int).tolist(),'zoom':c.astype(int),
             'psd_db':rounded(10*np.log10(np.maximum(ps[1:],1e-30)),3)}
    return features,aux,preview

def extract(args):
    args.cache.mkdir(parents=True,exist_ok=True)
    files=sorted(args.data.glob('*.cdf'))
    inventory=[]; rows=[]; features=[]; amplitudes=[]; previews=[]; countfeatures=[]
    if (args.cache/'features.npz').exists() and not args.refresh:
        old=np.load(args.cache/'features.npz')
        inventory=json.loads((args.cache/'inventory.json').read_text())
        if not {v['file'] for v in inventory}.issubset({f.name for f in files}):
            raise ValueError('Cached source file missing; investigate before refreshing')
        rows=json.loads((args.cache/'records.json').read_text())
        previews=json.loads((args.cache/'previews.json').read_text())
        features=[dict(zip(old['keys'].tolist(),v)) for v in old['X']]
        amplitudes=[dict(zip(old['akeys'].tolist(),v)) for v in old['A']]
        countfeatures=[dict(zip(old['ckeys'].tolist(),v)) for v in old['C']]
    previous={r['file']:r['sha256'] for r in inventory}
    meta=['Burst_ID','Burst_Type_Name','Burst_Quality','Burst_Quality_Algorithm_Name',
          'Trigger_Source_Name','Trigger_Position','Burst_Saturation_Flag','Sample_Speed',
          'N_Samples_per_Channel','Low_Pass_Filter','SC_Thrusting_Flag','SC_TWTA_On_Flag',
          'SCM_Cal_On_Flag','Burst_Time_Series_Physical_Units_Valid_Flag',
          'Burst_Total_SWEAP_Counts','SWEAP_Status','SWEAP_Electron_Mask','SWEAP_Ion_Mask','SWEAP_Start']
    for file in files:
        h=hashlib.sha256()
        with file.open('rb') as stream:
            for chunk in iter(lambda:stream.read(8*1024*1024), b''): h.update(chunk)
        if file.name in previous:
            if h.hexdigest()!=previous[file.name]:raise ValueError('Source changed; rerun with --refresh')
            print('Reusing verified extraction',file.name,flush=True)
            continue
        c=CDF(file); ep=c.varget('epoch'); nr=len(ep)
        metadata={v:c.varget(P+v) for v in meta}
        inventory.append({'file':file.name,'bytes':file.stat().st_size,'sha256':h.hexdigest(),
                          'records':nr,'variables':len(c.cdf_info().zVariables)})
        for rec in range(nr):
            row={'event':len(rows),'file':file.name,'record':rec,'tt2000':str(int(ep[rec])),
                 'utc':str(cdfepoch.to_datetime(ep[rec:rec+1])[0]),'day':file.name.split('_')[-2]}
            for key,vals in metadata.items():
                val=vals[rec]
                row[key]=str(val).strip() if isinstance(val,(str,np.str_)) else val.item()
            n=int(row['N_Samples_per_Channel']);fs=float(row['Sample_Speed'])
            names=channel_names(c,rec)
            row['analog_channels']=names
            row['analog_configuration']='reference' if names==CHANNELS else 'incomplete' if None in names else 'alternate_configuration'
            x=waveform(c,rec,n)
            fe,am,aux,preview=describe(x,fs)
            counts=c.varget(P+'Burst_Time_Series_SWEAP_Counts',startrec=rec,endrec=rec).ravel()[:n].astype(float)
            cf,ca,cp=describe_counts(counts,x,fs)
            if int(counts.sum())!=int(row['Burst_Total_SWEAP_Counts']):raise ValueError('Count total mismatch')
            cp['zoom']=cp['zoom'][preview['zoom_start']:preview['zoom_start']+1024].tolist()
            # Compact count spectrum on the same logarithmic bins as the analog PSDs.
            ff=np.arange(1,2049)*fs/4096;ee=np.unique(np.rint(np.geomspace(1,len(ff),180)).astype(int))-1
            cps=np.array(cp['psd_db']);cp['psd_db']=[float(np.mean(cps[l:r])) for l,r in zip(ee[:-1],ee[1:]) if r>l]
            preview['channels']=names
            if names!=CHANNELS:
                fe={k:float('nan') for k in fe};am={k:float('nan') for k in am}
                aux={k:None for k in aux}
                for k in cf:
                    if k.startswith('count_pair_'):cf[k]=float('nan')
                ca['counts_envelope_max_abs_corr']=None;ca['counts_shift_p_exploratory']=None
            preview['counts']=cp;row.update(ca);countfeatures.append(cf)
            row.update(aux);rows.append(row);features.append(fe);amplitudes.append(am);previews.append(preview)
            if rec%100==0:print(file.name,rec,'/',nr,flush=True)
    keys=list(features[0]);akeys=list(amplitudes[0])
    X=np.array([[r[k] for k in keys] for r in features]); A=np.array([[r[k] for k in akeys] for r in amplitudes])
    reference=np.array([r.get('analog_configuration','reference')=='reference' for r in rows])
    if not np.all(np.isfinite(X[reference])):raise ValueError('Nonfinite reference features')
    dump(args.cache/'records.json',rows);dump(args.cache/'previews.json',previews)
    dump(args.cache/'inventory.json',inventory)
    ckeys=list(countfeatures[0]);C=np.array([[r[k] for k in ckeys] for r in countfeatures])
    np.savez(args.cache/'features.npz',X=X,A=A,C=C,keys=keys,akeys=akeys,ckeys=ckeys)
    print('Extracted',len(rows),'bursts;',X.shape[1],'shape features',flush=True)

def robust_fit(X):
    mid=np.median(X,axis=0); scale=np.subtract(*np.percentile(X,[75,25],axis=0))
    scale=np.where(scale>1e-8,scale,np.std(X,axis=0))
    return mid,np.maximum(scale,1e-8)

def transform(X,mid,scale,keys):
    X=np.nan_to_num(np.clip((X-mid)/scale,-5,5),nan=0.) # excluded records: display projection only
    # Equal total squared-distance weight for temporal, spectral, coupling families.
    groups=[]
    for k in keys:
        if k.startswith('pair_'):g='coupling'
        elif any(s in k for s in ['spectral','log_f50','log_bandwidth','zero_crossing','autocorr']):g='spectral'
        elif 'log_rms' in k:g='amplitude'
        else:g='temporal'
        groups.append(g)
    weights=np.array([1/np.sqrt(groups.count(g)) for g in groups])
    return X*weights,groups

def labels_for(X,k,method='ward'):
    return cut_tree(linkage(X,method=method),n_clusters=k).ravel()

def silhouettes(D,labels):
    out=np.zeros(len(labels))
    for g in np.unique(labels):
        ii=np.flatnonzero(labels==g)
        if len(ii)<2:continue
        a=D[np.ix_(ii,ii)].sum(axis=1)/(len(ii)-1)
        b=np.min([D[np.ix_(ii,np.flatnonzero(labels==h))].mean(axis=1)
                  for h in np.unique(labels) if h!=g],axis=0)
        out[ii]=(b-a)/np.maximum(np.maximum(a,b),1e-30)
    return out

def ari(a,b):
    _,ia=np.unique(a,return_inverse=True);_,ib=np.unique(b,return_inverse=True)
    tab=np.zeros((ia.max()+1,ib.max()+1),int);np.add.at(tab,(ia,ib),1)
    choose=lambda v:np.sum(v*(v-1)/2)
    s=choose(tab);aa=choose(tab.sum(axis=1));bb=choose(tab.sum(axis=0));total=len(a)*(len(a)-1)/2
    expected=aa*bb/max(total,1);den=(aa+bb)/2-expected
    return float((s-expected)/den) if den else 1.

def align(a,b,k):
    tab=np.zeros((k,k),int);np.add.at(tab,(a,b),1)
    rr,cc=linear_sum_assignment(-tab);mapping=dict(zip(cc,rr))
    return np.array([mapping[v] for v in b])

def cluster(args):
    z=np.load(args.cache/'features.npz');X=z['X'];A=z['A'];keys=z['keys'].tolist()
    rows=json.loads((args.cache/'records.json').read_text())
    valid=np.array([not r['Burst_Saturation_Flag'] and r.get('analog_configuration','reference')=='reference' for r in rows])
    ix=np.flatnonzero(valid); mid,scale=robust_fit(X[ix]);Z,families=transform(X,mid,scale,keys)
    D=squareform(pdist(Z[ix]));tree=linkage(Z[ix],method='ward')
    candidates=[]
    for k in range(2,11):
        lab=cut_tree(tree,n_clusters=k).ravel();si=silhouettes(D,lab)
        candidates.append({'k':k,'silhouette':float(si.mean()),'sizes':np.bincount(lab).tolist()})
    eligible=[r for r in candidates if min(r['sizes'])>=20]
    best=max(eligible,key=lambda r:r['silhouette']);k=best['k']
    base=cut_tree(tree,n_clusters=k).ravel();sil=silhouettes(D,base)
    print('Candidate solutions:',json.dumps(candidates),flush=True)
    print('Selected',k,flush=True)
    rng=np.random.default_rng(SEED);agreement=np.zeros(len(ix));seen=np.zeros(len(ix));stability=[]
    for rep in range(50):
        sub=np.sort(rng.choice(len(ix),int(.8*len(ix)),replace=False))
        # Perturb both event sample and available feature set; refit scaling.
        cols=np.sort(rng.choice(X.shape[1],int(.8*X.shape[1]),replace=False))
        mm,ss=robust_fit(X[ix[sub]][:,cols]);zz,_=transform(X[ix[sub]][:,cols],mm,ss,[keys[j] for j in cols])
        lab=labels_for(zz,k);mapped=align(base[sub],lab,k)
        agreement[sub]+=(mapped==base[sub]);seen[sub]+=1;stability.append(ari(base[sub],lab))
    agreement/=np.maximum(seen,1)
    kth=np.sort(D,axis=1)[:,15];outlier=kth>np.quantile(kth,.98)
    # Avoid forcing tenuous assignments; thresholds are inspection policies, not probabilities.
    status=np.array(['core']*len(ix),dtype=object)
    status[(sil<.05)|(agreement<.8)]='ambiguous'
    status[outlier]='outlier'
    alts={}
    for title,cols in {
        'without_SCM':[j for j,key in enumerate(keys) if 'SCM5' not in key],
        'without_frequency_location':[j for j,key in enumerate(keys) if not any(w in key for w in ['log_f50','log_bandwidth','zero_crossing'])],
        'temporal_only':[j for j,g in enumerate(families) if g=='temporal'],
        'spectral_only':[j for j,g in enumerate(families) if g=='spectral'],
    }.items():
        zz,_=transform(X[ix][:,cols],mid[cols],scale[cols],[keys[j] for j in cols])
        alt=labels_for(zz,k);alts[title]={'ARI':ari(base,alt)}
    am,ass=robust_fit(A[ix]);za=np.clip((A-am)/ass,-5,5)/np.sqrt(A.shape[1])
    amp=labels_for(np.c_[Z[ix],za[ix]],k);alts['with_absolute_amplitude']={'ARI':ari(base,amp)}
    avg=labels_for(Z[ix],k,'average');alts['average_linkage']={'ARI':ari(base,avg),'sizes':np.bincount(avg).tolist()}
    # Independent per-day fits, compared only against that day's pooled assignments.
    for day in sorted(set(r['day'] for r in rows)):
        sub=np.flatnonzero(np.array([rows[j]['day']==day for j in ix]))
        if len(sub)<max(k,2):continue
        mm,ss=robust_fit(X[ix[sub]]);zz,_=transform(X[ix[sub]],mm,ss,keys)
        alt=labels_for(zz,k);alts['fit_day_'+day]={'ARI':ari(base[sub],alt),'n':len(sub)}
    centers=np.array([Z[ix][base==g].mean(axis=0) for g in range(k)])
    all_labels=np.argmin(cdist(Z,centers),axis=1);all_labels[ix]=base
    # Flagged records receive a nearest-centroid suggestion only; never a fitted class.
    for j,r in enumerate(rows):
        r['group']=int(all_labels[j]);r['status']='saturated';r['silhouette']=None;r['stability']=None
    for r in rows:
        if r.get('analog_configuration','reference')!='reference':r.update(group=-1,status=r['analog_configuration'])
    for a,j in enumerate(ix):
        rows[j].update(status=status[a],silhouette=float(sil[a]),stability=float(agreement[a]),
                       neighbor_distance=float(kth[a]))
    _,sv,V=np.linalg.svd(Z[ix]-Z[ix].mean(axis=0),full_matrices=False)
    pc=(Z-Z[ix].mean(axis=0))@V[:2].T
    for j,r in enumerate(rows):r.update(pc1=float(pc[j,0]),pc2=float(pc[j,1]))
    summaries=[]
    for g in range(k):
        ids=ix[base==g];dist=cdist(Z[ids],Z[ids]).sum(axis=1)
        ordered=ids[np.argsort(dist)]
        medoid=int(ordered[0]); core=[int(i) for i in ordered if rows[i]['status']=='core']
        # Include central examples, a far member, and the weakest assignment.
        representatives=list(dict.fromkeys([medoid]+core[:3]+[int(ids[np.argmax(cdist(Z[ids],centers[g:g+1]).ravel())]),int(ids[np.argmin(sil[base==g])])]))
        summaries.append({'group':g,'n':len(ids),'core':sum(rows[i]['status']=='core' for i in ids),
                          'medoid':medoid,'representatives':representatives,
                          'median_silhouette':float(np.median(sil[base==g])),
                          'median_stability':float(np.median(agreement[base==g])),
                          'feature_medians':dict(zip(keys,np.median(X[ids],axis=0).tolist()))})
    results={'method':'Ward hierarchical clustering in robust-scaled, equally weighted feature families',
             'seed':SEED,'k':k,'candidates':candidates,'silhouette':float(sil.mean()),
             'stability_ARI':{'median':float(np.median(stability)),'p10':float(np.quantile(stability,.1)),
                              'p90':float(np.quantile(stability,.9)),'values':stability},
             'sensitivity':alts,'groups':summaries,'pca_explained':(sv[:2]**2/np.sum(sv**2)).tolist(),
             'feature_names':keys,'feature_families':families,'scaling_median':mid.tolist(),'scaling_iqr':scale.tolist(),
             'centroids':centers.tolist(),'fit_indices':ix.tolist(),
             'outlier_threshold':float(np.quantile(kth,.98)),
             'policy':{'core_min_silhouette':.05,'core_min_resample_agreement':.8,'outlier_quantile':.98,'neighbor_k':15}}
    args.out.mkdir(parents=True,exist_ok=True)
    dump(args.cache/'classified.json',rows);dump(args.out/'model_diagnostics.json',results)
    with (args.out/'event_catalog.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));wr.writeheader();wr.writerows(rows)
    with (args.out/'features.csv').open('w',newline='') as f:
        wr=csv.writer(f);wr.writerow(['event']+keys+z['akeys'].tolist())
        for i in range(len(X)):wr.writerow([i]+X[i].tolist()+A[i].tolist())
    dump(args.out/'source_inventory.json',json.loads((args.cache/'inventory.json').read_text()))
    print('Stability and sensitivity:',json.dumps({'stability':results['stability_ARI']['median'],'sensitivity':alts}),flush=True)
    print('Status:',dict(zip(*np.unique([r['status'] for r in rows],return_counts=True))),flush=True)
    for g in summaries:
        med=g['feature_medians'];print('GROUP',g['group'],'n',g['n'],'core',g['core'],'medoid',g['medoid'],
          {ch:{a:round(med[ch+'__'+a],3) for a in ['log_crest','energy_width','spectral_entropy','log_f50','spectral_peak3']} for ch in CHANNELS},flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True)
    ap.add_argument('--cache',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--stage',choices=['extract','cluster','all'],default='all');ap.add_argument('--refresh',action='store_true');args=ap.parse_args()
    if args.stage in ['extract','all']:extract(args)
    if args.stage in ['cluster','all']:cluster(args)
