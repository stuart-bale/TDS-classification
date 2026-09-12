"""Exploratory count/voltage correlations, preserving burst and scan structure."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy import stats,signal,ndimage
from cdflib import CDF
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_tds import waveform,P,dump,CHANNELS
from build_atlas import savecsv,COLORS

def norm(a):
    a=a-a.mean(axis=-1,keepdims=True)
    return a/np.maximum(np.linalg.norm(a,axis=-1,keepdims=True),1e-30)

def linear_lags(a,b,lags):
    # Positive lag: counts b follow voltage a. Recenter each overlapping segment.
    result=[]
    for lag in lags:
        aa=a[:,:-lag] if lag>0 else a[:,-lag:] if lag<0 else a
        bb=b[lag:] if lag>0 else b[:lag] if lag<0 else b
        result.append(norm(aa)@norm(bb))
    return np.array(result).T

def run(args):
    rows=json.loads((args.cache/'final_records.json').read_text());prev=json.loads((args.cache/'previews.json').read_text())
    eligible=np.array([r.get('analog_configuration','reference')=='reference' and not r['Burst_Saturation_Flag'] for r in rows])
    ix=np.flatnonzero(eligible); within=np.array([i for i in ix if rows[i]['counts_total']>=32 and rows[i].get('counts_phase_baseline_status_matched',True)])
    # Cache 128-sample voltage means and RMS, with immutable source/record keys.
    ids=[f"{r['file']}:{r['record']}" for r in rows]
    path=args.cache/'voltage_bins.npz'
    if path.exists() and np.load(path)['ids'].tolist()==ids and 'native' in np.load(path):
        a=np.load(path);means=a['means'];env=a['env'];native=a['native']
    else:
        means=np.zeros((len(rows),4,256));env=np.zeros_like(means);native=np.full((len(rows),4,3),np.nan);last=None;cached=0
        if path.exists():
            old=np.load(path);oldids=old['ids'].tolist()
            if 'native' in old and oldids==ids[:len(oldids)]:
                cached=len(oldids);means[:cached]=old['means'];env[:cached]=old['env'];native[:cached]=old['native']
        for j in ix[ix>=cached]:
            r=rows[j]
            if r['file']!=last:c=CDF(args.data/r['file']);last=r['file']
            x=waveform(c,r['record'],r['N_Samples_per_Channel'])[:4]
            x-=np.median(x,axis=1,keepdims=True)
            if r['counts_total']>=32:
                cnt=c.varget(P+'Burst_Time_Series_SWEAP_Counts',startrec=r['record'],endrec=r['record']).ravel()[:r['N_Samples_per_Channel']].astype(float)
                native[j,:,0]=norm(x)@norm(cnt)
                freq,px=signal.welch(x,r['Sample_Speed'],nperseg=4096,noverlap=2048,axis=1)
                _,pc=signal.welch(cnt,r['Sample_Speed'],nperseg=4096,noverlap=2048)
                _,cross=signal.csd(x,cnt[None,:],r['Sample_Speed'],nperseg=4096,noverlap=2048,axis=1)
                coherence=abs(cross)**2/np.maximum(px*pc,1e-30)
                band=(freq>=1000)&(freq<=864000);hit=np.argmax(coherence[:,band],axis=1)
                native[j,:,1]=np.max(coherence[:,band],axis=1);native[j,:,2]=freq[band][hit]
            x=x.reshape(4,256,128)
            means[j]=x.mean(axis=2);env[j]=np.sqrt(np.mean(x*x,axis=2))
        np.savez(path,ids=ids,means=means,env=env,native=native)
    lags=np.arange(-16,17);fs=rows[0]['Sample_Speed'];lagms=lags*128/fs*1000
    rawcurves=[];rescurves=[];catalog=[];countbins=np.array([p['counts']['bins'] for p in prev],float)
    residual=np.array([p['counts']['phase_residual'] for p in prev],float)
    for j in within:
        # RMS captures modulation of high-frequency voltages; signed means only
        # describe a boxcar-filtered, coarse low-frequency voltage trace.
        a=np.r_[env[j],means[j]];b=countbins[j];rb=residual[j]
        raw=linear_lags(a,b,lags);res=linear_lags(a,rb,lags)
        rawcurves.append(raw);rescurves.append(res)
        typ,ch,li=0,0,0
        hit=np.unravel_index(np.argmax(abs(res)),res.shape);channel=hit[0]%4;kind='envelope' if hit[0]<4 else 'signed_mean'
        # Circular-shift reference: recompute max over all 8 traces and 33 lags.
        # This is a stationarity diagnostic, not a calibrated detection p-value.
        circ=np.fft.irfft(np.conj(np.fft.rfft(norm(a),axis=1))*np.fft.rfft(norm(rb)),n=256,axis=1)
        local=np.array([np.max(abs(circ[:,(lags+s)%256])) for s in range(256)])
        shift_p=float((1+np.sum(local[33:224]>=local[0]))/192)
        entry={k:rows[j][k] for k in ['event','file','record','utc','day','counts_total','SWEAP_Start','solar_distance_Rs']}
        entry.update(max_residual_abs_r=float(abs(res[hit])),residual_r_at_peak=float(res[hit]),peak_channel=CHANNELS[channel],peak_measure=kind,counts_lag_ms=float(lagms[hit[1]]),circular_shift_p_exploratory=shift_p)
        da=signal.detrend(a,axis=1);db=signal.detrend(rb)
        dr=linear_lags(da,db,lags);diff=norm(np.diff(a,axis=1))@norm(np.diff(rb))
        entry['peak_r_after_linear_detrend']=float(dr[hit])
        entry['max_abs_r_after_linear_detrend']=float(np.max(abs(dr)))
        entry['peak_trace_first_difference_r0']=float(diff[hit[0]])
        for c,ch in enumerate(CHANNELS[:4]):
            entry[f'{ch}_native_sample_r0']=float(native[j,c,0])
            entry[f'{ch}_max_welch_coherence']=float(native[j,c,1])
            entry[f'{ch}_peak_coherence_hz']=float(native[j,c,2])
            for k,offset in [('envelope',0),('signed_mean',4)]:
                entry[f'{ch}_{k}_raw_r0']=float(raw[c+offset,16]);entry[f'{ch}_{k}_residual_r0']=float(res[c+offset,16])
        catalog.append(entry)
    rawcurves=np.array(rawcurves);rescurves=np.array(rescurves)
    # Within-day, similar-phase and similar-total donor profiles preserve more
    # scan structure than arbitrary shifts. Leave target event out. This reports
    # reference percentiles only; sparse matching is explicitly unavailable.
    for q,j in enumerate(within):
        r=rows[j];donors=np.array([i for i in within if i!=j and rows[i]['day']==r['day'] and rows[i]['SWEAP_Status']==r['SWEAP_Status'] and abs(rows[i]['SWEAP_Start']-r['SWEAP_Start'])<=.01 and .25<=rows[i]['counts_total']/r['counts_total']<=4],int)
        catalog[q]['matched_donors']=len(donors)
        catalog[q]['matched_zero_lag_percentile']=None
        if len(donors)>=20:
            a=norm(np.r_[env[j],means[j]]);donorstat=np.max(abs(a@norm(residual[donors]).T),axis=0)
            observed=float(np.max(abs(rescurves[q,:,16])))
            catalog[q]['matched_zero_lag_percentile']=float(np.mean(donorstat<=observed)*100)
    # Across-burst ranks: amplitudes in engineering units, count totals, controls
    # withheld from the original clustering. No individual sample p-values.
    total=np.array([rows[i]['counts_total'] for i in ix]);phase=np.array([rows[i]['SWEAP_Start'] for i in ix]);radius=np.array([rows[i]['solar_distance_Rs'] for i in ix]);days=np.array([rows[i]['day'] for i in ix])
    t=(phase-phase.mean())/max(phase.std(),1e-20)
    design=[np.ones(len(ix)),t,t*t,t*t*t]
    for knot in np.quantile(t,[.1,.25,.5,.75,.9]):design.append(np.maximum(t-knot,0)**3)
    design.append(stats.rankdata(radius)/len(ix))
    for field in ['day','SWEAP_Status','Burst_Type_Name','SC_TWTA_On_Flag']:
        vals=np.array([str(rows[i][field]) for i in ix])
        design.extend((vals==v).astype(float) for v in np.unique(vals)[1:])
    design=np.array(design).T
    def remove(y):return y-design@np.linalg.lstsq(design,y,rcond=None)[0]
    yr=remove(stats.rankdata(total));across=[]
    for ch in CHANNELS[:4]:
        amp=np.array([rows[i][ch+'__rms'] for i in ix]);xr=remove(stats.rankdata(amp))
        across.append({'channel':ch,'n':len(ix),'pooled_spearman':float(stats.spearmanr(amp,total).statistic),'adjusted_rank_correlation':float(np.corrcoef(xr,yr)[0,1]),'by_day':{d:{'n':int(np.sum(days==d)),'spearman':float(stats.spearmanr(amp[days==d],total[days==d]).statistic)} for d in np.unique(days)}})
    summaries=[]
    for j in range(8):
        summaries.append({'channel':CHANNELS[j%4],'measure':'envelope' if j<4 else 'signed_mean','raw_zero_lag_median':float(np.median(rawcurves[:,j,16])),'residual_zero_lag_median':float(np.median(rescurves[:,j,16])),'residual_zero_lag_abs_ge_0_5':int(np.sum(abs(rescurves[:,j,16])>=.5)),'raw_median_lag_curve':np.median(rawcurves[:,j,:],axis=0).tolist(),'residual_median_lag_curve':np.median(rescurves[:,j,:],axis=0).tolist()})
    selected=sorted([r for r in catalog if r['counts_total']>=128],key=lambda r:r['max_residual_abs_r'],reverse=True)[:8]
    diag={'eligible_analog_bursts':len(ix),'within_burst_minimum_counts':32,'within_burst_n':len(within),'bin_samples':128,'bin_us':128/fs*1e6,'counts_positive_lag_definition':'Counts follow voltage','lags_ms':lagms.tolist(),'across_burst':across,'adjustment':'Rank residual correlation after day, SWEAP status, TWTA flag, selection type, radius rank and cubic spline of sweep phase (five interior knots). Descriptive observational adjustment, not a causal model.','within_burst':summaries,'matched_donor_eligible_n':sum(r['matched_donors']>=20 for r in catalog),'candidate_events':selected,'caveats':['SWEAP counts are scanned; phase residuals can retain scan structure.','Correlation coefficients quantify association, not wave-particle coupling.','Lag-peak selection searches eight traces and 33 lags; peak r is upward selected.','Circular shifts assume stationarity and matched donor profiles are exchangeable only approximately; neither is a calibrated physical-detection test.','Signed voltage uses 128-sample boxcar means, not resolved high-frequency voltage phase.','No timing offset or transfer-function correction is applied; lag accuracy is limited by binning and instrument timing.']}
    diag['native_cadence']={'sample_rate':fs,'n':len(within),'coherence_band_hz':[1000,864000],'welch_nperseg':4096,'welch_noverlap':2048,'channels':[{'channel':ch,'median_sample_r0':float(np.median(native[within,c,0])),'max_abs_sample_r0':float(np.max(abs(native[within,c,0]))),'median_peak_coherence':float(np.median(native[within,c,1]))} for c,ch in enumerate(CHANNELS[:4])],'caveat':'Native-sample correlations and peak coherence are uncorrected engineering diagnostics. Maxima search many frequency bins with only 15 overlapping Welch segments; no detection threshold or corrected phase delay is inferred.'}
    dump(args.out/'electron_voltage_diagnostics.json',diag);savecsv(args.out/'electron_voltage_correlations.csv',catalog)
    fig,axs=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    x=np.arange(4);axs[0,0].bar(x-.18,[a['pooled_spearman'] for a in across],width=.36,label='Pooled');axs[0,0].bar(x+.18,[a['adjusted_rank_correlation'] for a in across],width=.36,label='Adjusted');axs[0,0].set_xticks(x,CHANNELS[:4]);axs[0,0].set(ylabel='Rank correlation',title='Burst total counts vs voltage RMS');axs[0,0].legend()
    for j,ch in enumerate(CHANNELS[:4]):
        axs[0,1].plot(lagms,np.median(rawcurves[:,j,:],axis=0),ls='--',color=COLORS[j],alpha=.7)
        axs[0,1].plot(lagms,np.median(rescurves[:,j,:],axis=0),color=COLORS[j],label=ch)
        axs[1,0].hist(rescurves[:,j,16],bins=np.linspace(-1,1,41),histtype='step',color=COLORS[j],label=ch)
    axs[0,1].set(xlabel='Counts lag voltage (ms)',ylabel='Median within-burst r',title='Voltage envelope: raw dashed / phase residual solid');axs[0,1].legend()
    axs[1,0].set(xlabel='Zero-lag count-residual / voltage-envelope r',ylabel='Bursts',title='Individual bursts can differ from the median');axs[1,0].legend()
    table=np.array([[a['by_day'][d]['spearman'] for d in sorted(set(days))] for a in across]);im=axs[1,1].imshow(table,vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto');axs[1,1].set_xticks(range(len(set(days))),[d[-2:] for d in sorted(set(days))]);axs[1,1].set_yticks(range(4),CHANNELS[:4]);axs[1,1].set_title('Across-burst Spearman within each day');axs[1,1].set_xlabel('March 2025')
    for i in range(4):
        for j in range(table.shape[1]):axs[1,1].text(j,i,f'{table[i,j]:.2f}',ha='center',va='center')
    fig.colorbar(im,ax=axs[1,1],label='Spearman r')
    for ax in axs.ravel():ax.axhline(0,color='#999',lw=.5) if ax is not axs[1,1] else None
    fig.suptitle(f'SPAN-electron counts and voltage waveforms · {len(rows):,} recorded bursts\n{len(ix):,} common-configuration unsaturated bursts; {len(within):,} with ≥32 counts for within-burst analysis',fontsize=13)
    fig.savefig(args.figures/'electron_voltage_correlations.png',dpi=160);plt.close(fig)
    for item in selected[:4]:
        j=item['event'];q=np.flatnonzero(within==j)[0];fig,axs=plt.subplots(3,1,figsize=(12,8),layout='constrained');tms=(np.arange(256)+.5)*128/fs*1000
        axs[0].step(tms,countbins[j],where='mid',label='Observed counts');axs[0].plot(tms,prev[j]['counts']['phase_expected'],label='Phase baseline × event gain');axs[0].set_ylabel('Counts / 66.7 µs');axs[0].legend()
        for c,ch in enumerate(CHANNELS[:4]):axs[1].plot(tms,norm(env[j])[c]*16,label=ch,color=COLORS[c])
        selected_ch=CHANNELS.index(item['peak_channel'])
        if item['peak_measure']=='signed_mean':axs[1].plot(tms,norm(means[j])[selected_ch]*16,color='#666',ls='--',lw=1.5,label=item['peak_channel']+' signed mean')
        axs[1].plot(tms,norm(residual[j])*16,color='black',lw=1.5,label='Count residual');axs[1].set_ylabel('Centered / standard deviation');axs[1].legend(ncol=5)
        for c,ch in enumerate(CHANNELS[:4]):axs[2].plot(lagms,rescurves[q,c],color=COLORS[c],label=ch)
        if item['peak_measure']=='signed_mean':axs[2].plot(lagms,rescurves[q,4+selected_ch],color='#666',ls='--',lw=1.5,label=item['peak_channel']+' signed mean')
        axs[2].set(xlabel='Counts lag voltage (ms)',ylabel='Residual / voltage-trace r');axs[2].legend(ncol=4)
        fig.suptitle(f"Event {j} · {rows[j]['utc']} UTC\nSelected by largest residual correlation across voltage envelopes / signed means and lags; origin unassigned",fontsize=11)
        fig.savefig(args.figures/f'electron_voltage_event_{j:04d}.png',dpi=150);plt.close(fig)
    lines=['# Electron counts and voltage associations','',f'All {len(set(r["file"] for r in rows))} files: {len(rows):,} bursts. Across-burst comparisons use {len(ix):,} unsaturated records with the common analog configuration. Within-burst comparisons require at least 32 counts and a matching-status phase baseline ({len(within):,} bursts).','','| Voltage | Pooled count-total / RMS Spearman | Adjusted rank correlation |','|---|---:|---:|']
    lines += [f"| {a['channel']} | {a['pooled_spearman']:.3f} | {a['adjusted_rank_correlation']:.3f} |" for a in across]
    lines+=['',diag['adjustment'],'','| Voltage measure | Median raw zero-lag r | Median phase-residual zero-lag r | Number with residual absolute r ≥ 0.5 |','|---|---:|---:|---:|']
    lines += [f"| {s['channel']} {s['measure']} | {s['raw_zero_lag_median']:.3f} | {s['residual_zero_lag_median']:.3f} | {s['residual_zero_lag_abs_ge_0_5']} |" for s in summaries]
    lines+=['','Positive lag means counts follow voltage. The lag search covers ±1.067 ms in 66.7 µs steps. Voltage envelopes are RMS in 128-sample bins; signed means average those same samples. The high-frequency waveform phase is not resolved by this coarse signed trace.','','Candidates below have ≥128 counts and are ranked by their largest absolute phase-residual correlation across four envelopes, four signed means and 33 lags. This selection inflates peak correlations; these are inspection candidates, not detections.','','| Event | Counts | Channel / measure | Peak r | Counts lag (ms) | Matched donors |','|---|---:|---|---:|---:|---:|']
    lines += [f"| {r['event']} | {r['counts_total']} | {r['peak_channel']} / {r['peak_measure']} | {r['residual_r_at_peak']:.3f} | {r['counts_lag_ms']:.3f} | {r['matched_donors']} |" for r in selected]
    lines+=['',f"Matched donor comparisons were available for {diag['matched_donor_eligible_n']} events (at least 20 other bursts on the same day and SWEAP status, within 10 ms of sweep phase, with totals within a factor of four). The catalog retains zero-lag reference percentiles and exploratory circular-shift values. Neither establishes statistical significance under the instrument scan process.",'','The lag curves, distributions, daily comparisons and four candidate trace figures are saved in '+str(args.figures)+'. The full per-burst coefficients are in electron_voltage_correlations.csv.','']+[f'- {c}' for c in diag['caveats']]
    lines += ['','[Separate V5 / SCM4 analysis](alternate_configuration_report.md) covers the other analog configuration.','','## Slow-trend sensitivity','','At the originally selected trace and lag, removing a linear trend from both voltage and count residual gives:','','| Event | Original peak r | Detrended r at same trace/lag | First-difference zero-lag r, same trace |','|---|---:|---:|---:|']
    lines += [f"| {r['event']} | {r['residual_r_at_peak']:.3f} | {r['peak_r_after_linear_detrend']:.3f} | {r['peak_trace_first_difference_r0']:.3f} |" for r in selected]
    lines += ['','Differencing suppresses slow changes and amplifies counting noise; it probes timescale sensitivity, not whether the original association is physical.','','## Native-cadence check','','The catalog also contains direct signed voltage/count correlations at 1.92 MSa/s and peak magnitude-squared Welch coherence between 1 and 864 kHz. These use original engineering voltage samples and digital count increments. Spectra use 4096-sample windows with 50% overlap; frequency-response phase corrections and instrumental relative timing are unavailable. Peak coherence searches many bins and is not a detection statistic.','','| Voltage | Median native zero-lag r | Largest absolute native r | Median peak coherence |','|---|---:|---:|---:|']
    lines += [f"| {v['channel']} | {v['median_sample_r0']:.5f} | {v['max_abs_sample_r0']:.3f} | {v['median_peak_coherence']:.3f} |" for v in diag['native_cadence']['channels']]
    (args.out/'electron_voltage_report.md').write_text('\n'.join(lines))
    print(json.dumps({k:diag[k] for k in ['eligible_analog_bursts','within_burst_n','across_burst','matched_donor_eligible_n']}),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for name in ['data','cache','out','figures']:ap.add_argument('--'+name,type=Path,required=True)
    run(ap.parse_args())
