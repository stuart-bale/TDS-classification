"""Build six-channel atlas and count/joint diagnostics after analyze_tds.py.
Run with the same --data, --cache, --out arguments. See README.md.
"""
import argparse, base64, csv, html, json, sys
from pathlib import Path
from collections import Counter
import numpy as np
from scipy import stats, signal
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.cluster.hierarchy import linkage, cut_tree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cdflib import CDF
from analyze_tds import (P, CHANNELS, robust_fit, labels_for, silhouettes, ari, align,
                         dump, waveform, SEED, transform, channel_names)

COLORS=['#007c91','#d86b25','#7d4db0','#388a52','#bd4671','#a39520','#497dcc','#555555','#00a59b','#a0522d']

def savecsv(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def simple_fit(X, maxk=8, reps=40):
    mid,scale=robust_fit(X);Z=np.clip((X-mid)/scale,-5,5)/np.sqrt(X.shape[1]);D=squareform(pdist(Z))
    tree=linkage(Z,method='ward');candidates=[]
    for k in range(2,maxk+1):
        labels=cut_tree(tree,n_clusters=k).ravel();si=silhouettes(D,labels)
        candidates.append({'k':k,'silhouette':float(si.mean()),'sizes':np.bincount(labels).tolist()})
    best=max([r for r in candidates if min(r['sizes'])>=20],key=lambda r:r['silhouette'])
    k=best['k'];lab=cut_tree(tree,n_clusters=k).ravel();si=silhouettes(D,lab)
    agree=np.zeros(len(X));seen=np.zeros(len(X));scores=[];rng=np.random.default_rng(SEED+6)
    for rep in range(reps):
        idx=np.sort(rng.choice(len(X),int(.8*len(X)),replace=False))
        cols=np.sort(rng.choice(X.shape[1],max(3,int(.8*X.shape[1])),replace=False))
        mm,ss=robust_fit(X[idx][:,cols]);zz=np.clip((X[idx][:,cols]-mm)/ss,-5,5)
        ll=labels_for(zz,k);agree[idx]+=align(lab[idx],ll,k)==lab[idx];seen[idx]+=1;scores.append(ari(lab[idx],ll))
    agree/=np.maximum(seen,1)
    return lab,si,agree,Z,{'k':k,'candidates':candidates,'silhouette':float(si.mean()),
                         'resample_ARI_median':float(np.median(scores)),
                         'resample_ARI_p10':float(np.quantile(scores,.1)),
                         'feature_median':mid.tolist(),'feature_scale':scale.tolist()}

def counts_and_joint(args,rows,z,model):
    C=z['C'];X=z['X'];ckeys=z['ckeys'].tolist();keys=z['keys'].tolist()
    usable=np.array([r['counts_usable_shape'] for r in rows]);ix=np.flatnonzero(usable)
    nc=sum(k.startswith('counts__') for k in ckeys)
    lab,si,ag,zc,cd=simple_fit(C[ix,:nc]);cd['features']=ckeys[:nc];cd['minimum_counts']=32
    for r in rows:
        r['count_group']=-1;r['count_status']='channel absent' if not r.get('counts_present',True) else 'zero counts' if r['counts_total']==0 else 'too few counts'
        r['count_stability']=None;r['count_silhouette']=None
    groups=[]
    for j,i in enumerate(ix):
        rows[i].update(count_group=int(lab[j]),count_status='core' if si[j]>=.05 and ag[j]>=.8 else 'ambiguous',
                       count_stability=float(ag[j]),count_silhouette=float(si[j]))
    for g in range(cd['k']):
        loc=np.flatnonzero(lab==g);ids=ix[loc];med=int(ids[np.argmin(cdist(zc[loc],zc[loc]).sum(axis=1))])
        groups.append({'group':g,'n':len(ids),'medoid':med,'median_total':float(np.median([rows[i]['counts_total'] for i in ids])),
                       'median_sweep_phase_s':float(np.median([rows[i]['SWEAP_Start'] for i in ids])),
                       'core':sum(rows[i]['count_status']=='core' for i in ids),
                       'feature_medians':dict(zip(ckeys[:nc],np.median(C[ids,:nc],axis=0).tolist()))})
    cd['groups']=groups
    # Joint clustering: analog representation = one block, count shape = one block,
    # cross-count/analog envelope relationships = one block. Absolute count rate excluded.
    jj=np.array([i for i in ix if rows[i]['status'] in ['core','ambiguous','outlier']]);mid=np.array(model['scaling_median']);sc=np.array(model['scaling_iqr'])
    za,_=transform(X[jj],mid,sc,keys);cm,cs=robust_fit(C[jj,:nc]);cc=np.clip((C[jj,:nc]-cm)/cs,-5,5)/np.sqrt(nc)
    rm,rs=robust_fit(C[jj,nc:]);rel=np.clip((C[jj,nc:]-rm)/rs,-5,5)/np.sqrt(C.shape[1]-nc)
    joint=np.c_[za/np.sqrt(3),cc,rel]
    jk=model['k'];jl=labels_for(joint,jk);anal=np.array([rows[i]['group'] for i in jj]);jd={'k':jk,'n':len(jj),'ARI_vs_analog':ari(anal,jl)}
    jd['counts_vs_analog_ARI']=ari(anal,np.array([rows[i]['count_group'] for i in jj]))
    # No shared integer naming semantics: J0 is not analog G0.
    for r in rows:r['joint_group']=-1
    for j,i in enumerate(jj):rows[i]['joint_group']=int(jl[j])
    rate=np.log1p(np.array([r['counts_total'] for r in rows])[ix])[:,None];rm,rs=robust_fit(rate)
    withrate=labels_for(np.c_[zc,np.clip((rate-rm)/rs,-5,5)],cd['k'])
    cd['ARI_with_count_rate']=ari(lab,withrate)
    cd['status_counts']=dict(Counter(r['count_status'] for r in rows))
    # Diagnostics for instrumental/selection dependence, rather than using them as labels.
    nuisance={}
    for view,key in [('analog','group'),('counts','count_group')]:
        data=[r for r in rows if r[key]>=0 and (key!='group' or r['status'] in ['core','ambiguous','outlier'])]
        nuisance[view]={}
        for field in ['day','Burst_Type_Name','SWEAP_Status']:
            gs=sorted(set(r[key] for r in data));cats=sorted(set(str(r[field]) for r in data))
            tab=np.array([[sum(r[key]==g and str(r[field])==cat for r in data) for cat in cats] for g in gs])
            chi=stats.chi2_contingency(tab,correction=False)[0]
            v=np.sqrt(chi/(len(data)*max(1,min(tab.shape)-1)))
            nuisance[view][field]={'columns':cats,'rows':gs,'table':tab.tolist(),'cramers_v':float(v)}
    phase=np.array([r['SWEAP_Start'] for r in rows]);total=np.array([r['counts_total'] for r in rows])
    rho,p=stats.spearmanr(phase,total);cd['count_total_vs_sweep_phase_spearman']=float(rho)
    # BH q-values for the exploratory max-over-channel, max-over-lag shift screen.
    tested=np.flatnonzero(usable & np.array([r.get("analog_configuration","reference")=="reference" for r in rows]));pv=np.array([rows[i]['counts_shift_p_exploratory'] for i in tested]);order=np.argsort(pv)
    q=np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1];qq=np.empty(len(q));qq[order]=np.minimum(q,1)
    for r in rows:r['counts_shift_q_exploratory']=None
    for j,i in enumerate(tested):rows[i]['counts_shift_q_exploratory']=float(qq[j])
    cd['exploratory_shift_tests']={'n':len(tested),'BH_q_below_0_05':int(np.sum(qq<.05)),
       'caveat':'Circular-shift stationarity assumption can fail for scanned counts; association is not wave-particle coupling.'}
    dump(args.out/'counts_and_joint_diagnostics.json',{'counts':cd,'joint':jd,'nuisance':nuisance})
    savecsv(args.out/'event_catalog.csv',rows)
    savecsv(args.out/'count_features.csv',[dict(event=i,**dict(zip(ckeys,C[i].tolist()))) for i in range(len(rows))])
    dump(args.cache/'final_records.json',rows)
    print('COUNT GROUPS',json.dumps(cd),flush=True);print('JOINT',jd,flush=True);print('NUISANCE',json.dumps(nuisance),flush=True)
    return cd,jd,nuisance

def plot_event(args,r,filename):
    c=CDF(args.data/r['file']);n=r['N_Samples_per_Channel'];fs=r['Sample_Speed'];x=waveform(c,r['record'],n)
    x-=np.median(x,axis=1,keepdims=True);cnt=c.varget(P+'Burst_Time_Series_SWEAP_Counts',startrec=r['record'],endrec=r['record']).ravel()[:n]
    z=x/np.maximum(np.sqrt(np.mean(x*x,axis=1,keepdims=True)),1e-20);peak=int(np.argmax(np.mean(z[:4]**2,axis=0)))
    lo=int(np.clip(peak-512,0,n-1024));tt=np.arange(n)/fs*1000
    fig,axs=plt.subplots(6,3,figsize=(15,10),gridspec_kw={'width_ratios':[1.7,1.15,1.2]},layout='constrained')
    for j,ch in enumerate(channel_names(c,r['record'])):
        if ch is None:
            for ax in axs[j]:ax.text(.5,.5,'Channel absent',transform=ax.transAxes,ha='center')
            continue
        b=x[j].reshape(512,-1);t=tt.reshape(512,-1).mean(axis=1)
        axs[j,0].vlines(t,b.min(axis=1),b.max(axis=1),color=COLORS[j],lw=.65)
        axs[j,1].plot(tt[lo:lo+1024],x[j,lo:lo+1024],color=COLORS[j],lw=.65)
        f,ps=signal.welch(x[j],fs,nperseg=4096,noverlap=2048)
        axs[j,2].loglog(f[1:]/1000,ps[1:],color=COLORS[j],lw=.7)
        unit='nT' if j==4 else 'mV';axs[j,0].set_ylabel(ch+'\n'+unit)
        axs[j,2].set_ylabel(unit+'²/Hz',fontsize=8)
        axs[j,0].axvspan(tt[lo],tt[lo+1023],color='#888888',alpha=.12)
    bins=cnt.reshape(256,-1).sum(axis=1);t=tt.reshape(256,-1).mean(axis=1)
    axs[5,0].step(t,bins,where='mid',color='#222222',lw=.8);axs[5,0].set_ylabel('SPAN-e\ncounts/66.7 µs')
    axs[5,1].step(tt[lo:lo+1024],cnt[lo:lo+1024],where='mid',color='#222222',lw=.7)
    axs[5,1].set_ylabel('counts/sample',fontsize=8)
    f,ps=signal.welch(cnt,fs,nperseg=4096,noverlap=2048)
    if np.any(ps[1:]>0):axs[5,2].loglog(f[1:]/1000,ps[1:],color='#222222',lw=.7)
    else:axs[5,2].text(.5,.5,'Zero counts',transform=axs[5,2].transAxes,ha='center')
    axs[5,2].set_ylabel('counts²/Hz',fontsize=8)
    if not r.get('counts_present',True):
        for col in range(3):
            axs[5,col].clear();axs[5,col].text(.5,.5,'SPAN-e channel absent',transform=axs[5,col].transAxes,ha='center');axs[5,col].set_axis_off()
    for j in range(6):
        for col in range(3):axs[j,col].grid(alpha=.15);axs[j,col].tick_params(labelsize=8)
        axs[j,0].set_xlim(0,n/fs*1000)
    axs[0,0].set_title('Full record\nmin/max envelope',fontsize=10);axs[0,1].set_title('1024 actual samples\nat strongest electric excursion',fontsize=10);axs[0,2].set_title('Welch spectrum\n468.75 Hz bins',fontsize=10)
    axs[5,0].set_xlabel('Time since burst start (ms)');axs[5,1].set_xlabel('Time since burst start (ms)');axs[5,2].set_xlabel('Frequency (kHz)')
    label=f"G{r['group']} / {r['status']}" if r['group']>=0 else f"A{r['alternate_group']} / {r['alternate_status']}" if r.get('alternate_group',-1)>=0 else r['status']
    fig.suptitle(f"Event {r['event']} · {r['utc'][:26]} UTC · {label} · C{r['count_group']} / {r['count_status']}\n"
        f"Engineering units, median removed. SPAN-e total {format(r['counts_total'], ',') if r.get('counts_present',True) else 'unavailable'}; sweep phase {r['SWEAP_Start']:.6f} s; status {r['SWEAP_Status']}",fontsize=11)
    fig.savefig(filename,dpi=135);plt.close(fig)

def overview(args,rows,model,cd,nuisance):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,3,figsize=(15,9),layout='constrained')
    for g in range(model['k']):
        rr=[r for r in rows if r['group']==g and r['status'] in ['core','ambiguous','outlier']]
        axs[0,0].scatter([r['pc1'] for r in rr],[r['pc2'] for r in rr],s=12,alpha=.55,color=COLORS[g],label=f'G{g} ({len(rr)})')
    rr=[r for r in rows if r['status'] in ['saturated','outlier']]
    axs[0,0].scatter([r['pc1'] for r in rr],[r['pc2'] for r in rr],s=38,facecolor='none',edgecolor='#222',linewidth=.8)
    axs[0,0].set(xlabel=f"PC1 ({model['pca_explained'][0]:.0%} variance)",ylabel=f"PC2 ({model['pca_explained'][1]:.0%})",title='Analog morphology · projection only')
    axs[0,0].legend(fontsize=8)
    for g in range(model['k']):
        rr=[r for r in rows if r['group']==g and r['status'] in ['core','ambiguous','outlier']]
        axs[0,1].scatter([r['V3V4__f50_hz']/1000 for r in rr],[10**float(FEATURES[r['event'],KEYS.index('V3V4__log_crest')]) for r in rr],s=10,alpha=.5,color=COLORS[g])
    axs[0,1].set(xscale='log',yscale='log',xlabel='V3−V4 spectral median (kHz)',ylabel='V3−V4 peak / RMS',title='Frequency content and impulsiveness')
    cand=model['candidates'];axs[0,2].plot([r['k'] for r in cand],[r['silhouette'] for r in cand],'o-',color=COLORS[0])
    axs[0,2].axvline(model['k'],ls=':',color='#555');axs[0,2].set(xlabel='Number of groups',ylabel='Mean silhouette',title='Separation across candidate partitions')
    for g in sorted(set(r['count_group'] for r in rows)):
        rr=[r for r in rows if r['count_group']==g and r.get('counts_present',True)]
        axs[1,0].scatter([r['SWEAP_Start'] for r in rr],[r['counts_total']+1 for r in rr],s=11,alpha=.5,color=COLORS[g] if g>=0 else '#aaa',label=f'C{g}' if g>=0 else 'Insufficient counts')
    axs[1,0].set(yscale='log',xlabel='Time since SWEAP sweep start (s)',ylabel='Counts per burst + 1',title='Digital-channel groups and sweep phase');axs[1,0].legend(fontsize=8)
    tab=nuisance['analog']['Burst_Type_Name'];a=np.array(tab['table']);a=a/a.sum(axis=0,keepdims=True)
    im=axs[1,1].imshow(a,vmin=0,vmax=1,cmap='Blues',aspect='auto');axs[1,1].set_xticks(range(len(tab['columns'])),tab['columns']);axs[1,1].set_yticks(range(model['k']),[f'G{i}' for i in range(model['k'])]);axs[1,1].set_title('Analog groups within each selection type')
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):axs[1,1].text(j,i,f'{a[i,j]:.0%}',ha='center',va='center',color='white' if a[i,j]>.5 else 'black')
    st=Counter(r['status'] for r in rows);axs[1,2].bar(st.keys(),st.values(),color=['#007c91','#d86b25','#999','#bd4671'][:len(st)])
    axs[1,2].set(ylabel='Bursts',title='Assignment status');axs[1,2].set_xticks(range(len(st)),[k.replace('alternate_configuration','alternate').replace('incomplete','missing') for k in st]);axs[1,2].tick_params(axis='x',rotation=30)
    for i,(k,v) in enumerate(st.items()):axs[1,2].text(i,v+4,str(v),ha='center')
    fig.suptitle(f"Parker Solar Probe / TDS · {min(r['utc'] for r in rows)[:10]} to {max(r['utc'] for r in rows)[:10]} · {len(rows)} six-channel bursts",fontsize=16)
    fig.savefig(args.figures/'overview.png',dpi=160);plt.close(fig)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);ap.add_argument('--cache',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--figures',type=Path,required=True);ap.add_argument('--skip-counts',action='store_true');ap.add_argument('--skip-examples',action='store_true');args=ap.parse_args()
    args.figures.mkdir(parents=True,exist_ok=True)
    z=np.load(args.cache/'features.npz');FEATURES=z['X'];KEYS=z['keys'].tolist()
    model=json.loads((args.out/'model_diagnostics.json').read_text());rows=json.loads((args.cache/'classified.json').read_text())
    if args.skip_counts:
        rows=json.loads((args.cache/'final_records.json').read_text());d=json.loads((args.out/'counts_and_joint_diagnostics.json').read_text());cd,jd,nuisance=d['counts'],d['joint'],d['nuisance']
    else:cd,jd,nuisance=counts_and_joint(args,rows,z,model)
    overview(args,rows,model,cd,nuisance)
    selected=set()
    for g in model['groups']:selected.update(g['representatives'])
    selected.update(g['medoid'] for g in cd['groups'])
    selected.update(r['event'] for r in sorted([r for r in rows if r['status']=='outlier'],key=lambda r:r['neighbor_distance'],reverse=True)[:6])
    selected.update(r['event'] for r in rows if r['Burst_Saturation_Flag'] or r['status']=='incomplete')
    for day in sorted({r['day'] for r in rows}):selected.update(r['event'] for r in [v for v in rows if v['day']==day and v['status']=='alternate_configuration'][:2])
    (args.figures/'examples').mkdir(exist_ok=True)
    if not args.skip_examples:
        for i in sorted(selected):plot_event(args,rows[i],args.figures/'examples'/f'event_{i:04d}.png')
    dump(args.cache/'selected_examples.json',sorted(selected))
    print('Rendered',len(selected),'six-channel examples',flush=True)
