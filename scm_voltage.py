"""SCM/voltage association screens in engineering units, not physical mode labels."""
import argparse,json
from pathlib import Path
from collections import Counter
import numpy as np
from scipy import signal,ndimage
from cdflib import CDF
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_tds import waveform,channel_names,dump
from build_atlas import savecsv,COLORS
from electron_voltage import norm,linear_lags

def spectral(x,fs):
 f,ps=signal.welch(x,fs,nperseg=4096,noverlap=2048,axis=1)
 _,cross=signal.csd(x[:4],x[4:5],fs,nperseg=4096,noverlap=2048,axis=1)
 co=np.minimum(1,abs(cross)**2/np.maximum(ps[:4]*ps[4],1e-30))
 return f,ps,cross,co

def run(args):
 rows=json.loads((args.cache/'final_records.json').read_text());spec=np.load(args.cache/'spectra.npy',mmap_mode='r');fs=rows[0]['Sample_Speed'];freq=np.arange(1,2049)*fs/4096
 eligible=[i for i,r in enumerate(rows) if not r['Burst_Saturation_Flag'] and r.get('analog_configuration','reference')!='incomplete']
 states={'reference':['V2','V4','V1V2','V3V4','SCM5'],'alternate_configuration':['V2','V5','V1V2','V3V4','SCM4']};masks={};maskinfo={}
 for state,names in states.items():
  ix=[i for i in eligible if rows[i].get('analog_configuration','reference')==state];a=10*np.log10(np.maximum(np.asarray(spec[ix]),1e-30));excess=a-ndimage.median_filter(a,size=(1,1,31),mode='nearest');med=np.median(excess,axis=0);persistence=np.mean(excess>8,axis=0);mm=[];info=[]
  for ch,name in enumerate(names):
   peaks,_=signal.find_peaks(med[ch],height=8,distance=3);peaks=[v for v in peaks if persistence[ch,v]>=.6 and v>=16];mask=np.zeros(2049,bool)
   for v in peaks:mask[max(1,v):min(2049,v+3)]=True
   mm.append(mask);info.append({'channel':name,'line_centers_hz':freq[peaks].tolist()})
  masks[state]=np.array(mm);maskinfo[state]={'n_training_records':len(ix),'channels':info}
  del a,excess
 pairs=[];summaries=[];last=None;lags=np.arange(-16,17);lagms=lags*128/fs*1000
 for i,r in enumerate(rows):
  r['scm_voltage_flags']='not assessed';r['scm_voltage_max_outside_line_coherence']=None;r['scm_voltage_max_envelope_abs_r']=None
 for j in eligible:
  r=rows[j];state=r.get('analog_configuration','reference');names=states[state]
  if last!=r['file']:c=CDF(args.data/r['file']);last=r['file']
  x=waveform(c,r['record'],r['N_Samples_per_Channel']);x-=np.median(x,axis=1,keepdims=True);f,ps,cross,co=spectral(x,fs);band=(f>=1000)&(f<=864000)
  native=norm(x[:4])@norm(x[4]);env=np.sqrt(np.mean(x.reshape(5,256,128)**2,axis=2));ec=linear_lags(env[:4],env[4],lags);allflags=[];eventpairs=[]
  for ch,name in enumerate(names[:4]):
   mask=masks[state][ch]|masks[state][4];keep=band&~mask;weight=np.sqrt(ps[ch]*ps[4]);raw=float(np.sum(co[ch,band]*weight[band])/max(weight[band].sum(),1e-30));outside=float(np.sum(co[ch,keep]*weight[keep])/max(weight[keep].sum(),1e-30));line_weight=float(weight[band&mask].sum()/max(weight[band].sum(),1e-30));scm_remaining=float(ps[4,keep].sum()/max(ps[4,band].sum(),1e-30));peak=np.argmax(abs(ec[ch]));spectralpeak=np.flatnonzero(keep)[np.argmax(co[ch,keep])]
   runs=np.flatnonzero(np.diff(np.r_[False,keep&(co[ch]>=.6),False]));longest=int(max((b-a for a,b in zip(runs[::2],runs[1::2])),default=0))
   flags=[]
   if abs(native[ch])>=.5:flags.append('waveform correlation')
   if abs(ec[ch,peak])>=.5:flags.append('envelope association')
   if outside>=.3:flags.append('coherence outside persistent lines')
   if raw>=.3 and line_weight>=.8 and outside<.3:flags.append('persistent-line dominated coherence')
   item={k:r[k] for k in ['event','file','record','utc','day','SC_TWTA_On_Flag']};item.update(voltage=name,scm=names[4],native_r0=float(native[ch]),envelope_r0=float(ec[ch,16]),envelope_peak_r=float(ec[ch,peak]),scm_lag_ms=float(lagms[peak]),weighted_coherence_raw=raw,weighted_coherence_outside_lines=outside,masked_cross_spectral_weight_fraction=line_weight,scm_power_remaining_fraction=scm_remaining,peak_coherence_outside_lines=float(co[ch,spectralpeak]),peak_coherence_frequency_hz=float(f[spectralpeak]),longest_coherent_run_hz=longest*fs/4096,flags='; '.join(flags) if flags else 'no threshold crossed')
   pairs.append(item);eventpairs.append(item);allflags.extend(flags)
  rows[j]['scm_voltage_flags']='; '.join(sorted(set(allflags))) if allflags else 'no threshold crossed';rows[j]['scm_voltage_max_outside_line_coherence']=max(v['weighted_coherence_outside_lines'] for v in eventpairs);rows[j]['scm_voltage_max_envelope_abs_r']=max(abs(v['envelope_peak_r']) for v in eventpairs)
  if j%500==0:print('SCM associations',j,'/',len(rows),flush=True)
 savecsv(args.out/'scm_voltage_correlations.csv',pairs);dump(args.cache/'final_records.json',rows)
 for state,names in states.items():
  for name in names[:4]:
   rr=[r for r in pairs if r['voltage']==name and r['scm']==names[4]]
   summaries.append({'voltage':name,'scm':names[4],'n':len(rr),'median_native_r0':float(np.median([r['native_r0'] for r in rr])),'median_envelope_r0':float(np.median([r['envelope_r0'] for r in rr])),'median_raw_coherence':float(np.median([r['weighted_coherence_raw'] for r in rr])),'median_outside_line_coherence':float(np.median([r['weighted_coherence_outside_lines'] for r in rr])),'flag_counts':{flag:sum(flag in r['flags'] for r in rr) for flag in ['waveform correlation','envelope association','coherence outside persistent lines','persistent-line dominated coherence']}})
 byevent={}
 for r in sorted(pairs,key=lambda r:r['weighted_coherence_outside_lines'],reverse=True):byevent.setdefault(r['event'],r)
 selected=list(byevent.values())[:4]
 lineevents=sorted([r for r in pairs if 'persistent-line dominated coherence' in r['flags']],key=lambda r:r['weighted_coherence_raw'],reverse=True)
 if lineevents and lineevents[0]['event'] not in [r['event'] for r in selected]:selected.append(lineevents[0])
 diag={'n_records':len(rows),'n_assessed':len(eligible),'n_pairs':len(pairs),'thresholds':{'absolute_native_r0':.5,'absolute_peak_envelope_r':.5,'weighted_coherence_outside_lines':.3,'line_dominated_min_raw_coherence':.3,'line_dominated_min_cross_spectral_weight_fraction':.8},'frequency_band_hz':[1000,864000],'lag_ms':lagms.tolist(),'positive_lag':'SCM envelope follows voltage envelope','line_masks':maskinfo,'summaries':summaries,'candidate_events':selected,'event_flag_counts':{flag:sum(flag in r['scm_voltage_flags'] for r in rows) for flag in ['waveform correlation','envelope association','coherence outside persistent lines','persistent-line dominated coherence','no threshold crossed','not assessed']},'caveats':['Flags are descriptive screening thresholds, not significance tests or physical wave identifications.','SCM4 and SCM5 are evaluated separately; saturated or missing-SCM records are not assessed.','Persistent lines are learned separately per configuration (>8 dB local excess in >=60% of unsaturated records, +/-1 Welch bin); persistence does not establish instrumental origin.','Welch coherence uses 15 overlapping 4096-sample segments; individual spectral maxima are upward selected and are not independently significant.','Magnitude-squared coherence is averaged with sqrt(voltage PSD × SCM PSD) weights; removing lines changes both bins and weights.','No transfer-function phase correction or instrument-relative timing correction is applied. Coherence does not establish physical E/B polarization, propagation, impedance or Poynting flux.']}
 dump(args.out/'scm_voltage_diagnostics.json',diag)
 fig,axs=plt.subplots(2,2,figsize=(13,8),layout='constrained')
 for col,(state,names) in enumerate(states.items()):
  ss=[s for s in summaries if s['scm']==names[4]];xx=np.arange(4)
  axs[0,col].bar(xx-.17,[s['median_raw_coherence'] for s in ss],width=.34,label='All band bins');axs[0,col].bar(xx+.17,[s['median_outside_line_coherence'] for s in ss],width=.34,label='Outside persistent lines');axs[0,col].set_xticks(xx,names[:4]);axs[0,col].set(ylabel='Median weighted coherence',title=names[4]+' / voltage spectral association');axs[0,col].legend(fontsize=8)
  for ch,name in enumerate(names[:4]):axs[1,col].hist([r['envelope_r0'] for r in pairs if r['voltage']==name and r['scm']==names[4]],bins=np.linspace(-1,1,41),histtype='step',label=name,color=COLORS[ch])
  axs[1,col].set(xlabel='Zero-lag RMS-envelope correlation',ylabel='Bursts',title=names[4]+' / voltage envelopes');axs[1,col].legend(fontsize=8)
 fig.suptitle(f'SCM–voltage associations · {len(eligible):,} unsaturated complete records\nThreshold flags guide inspection; they do not identify physical modes',fontsize=13);fig.savefig(args.figures/'scm_voltage_associations.png',dpi=160);plt.close(fig)
 for item in selected:
  r=rows[item['event']];c=CDF(args.data/r['file']);x=waveform(c,r['record'],r['N_Samples_per_Channel']);x-=np.median(x,axis=1,keepdims=True);state=r.get('analog_configuration','reference');names=states[state];ch=names.index(item['voltage']);f,ps,cross,co=spectral(x,fs);env=np.sqrt(np.mean(x.reshape(5,256,128)**2,axis=2));ec=linear_lags(env[ch:ch+1],env[4],lags)[0];tt=(np.arange(256)+.5)*128/fs*1000
  fig,axs=plt.subplots(2,2,figsize=(13,8),layout='constrained');axs[0,0].plot(tt,norm(env[ch])*16,label=item['voltage']);axs[0,0].plot(tt,norm(env[4])*16,label=item['scm']);axs[0,0].set(xlabel='Time (ms)',ylabel='Envelope, centered / standard deviation');axs[0,0].legend()
  center=int(np.argmax(abs(x[ch])));lo=int(np.clip(center-512,0,len(x[ch])-1024));t=np.arange(lo,lo+1024)/fs*1000
  for idx,name in [(ch,item['voltage']),(4,item['scm'])]:axs[0,1].plot(t,x[idx,lo:lo+1024]/max(np.sqrt(np.mean(x[idx]**2)),1e-30),label=name,lw=.8)
  axs[0,1].set(xlabel='Time (ms) · 1024 actual samples',ylabel='Median removed / full-record RMS');axs[0,1].legend()
  band=(f>=1000)&(f<=864000);mask=masks[state][ch]|masks[state][4];axs[1,0].semilogx(f[band]/1000,co[ch,band],lw=.7,label='Coherence');axs[1,0].scatter(f[band&mask]/1000,co[ch,band&mask],s=9,color='#c44',label='Masked persistent-line bins');axs[1,0].set(xlabel='Frequency (kHz)',ylabel='Magnitude-squared coherence',ylim=(0,1));axs[1,0].legend(fontsize=8)
  axs[1,1].plot(lagms,ec);axs[1,1].set(xlabel='SCM envelope lags voltage (ms)',ylabel='Envelope correlation');fig.suptitle(f"Event {item['event']} · {r['utc']} UTC\n{item['voltage']} / {item['scm']} · {item['flags']}",fontsize=12);fig.savefig(args.figures/f"scm_voltage_event_{item['event']:04d}.png",dpi=150);plt.close(fig)
 lines=['# SCM–voltage association flags','',f"{len(eligible):,} of {len(rows):,} bursts assessed ({len(pairs):,} voltage/SCM pairs). Saturation and missing-SCM records remain unassessed. All four available voltages are tested against SCM5 or SCM4, keeping configurations separate.",'','Flags are deliberately descriptive and may overlap:','','- **Waveform correlation:** absolute full-sample zero-lag Pearson r ≥ 0.5.','- **Envelope association:** absolute peak RMS-envelope r ≥ 0.5 across ±1.067 ms (66.7 µs bins); positive lag means SCM follows voltage.','- **Coherence outside persistent lines:** PSD-weighted mean magnitude-squared coherence ≥ 0.3 after excluding learned narrow-line bins. This does not necessarily mean broadband coherence.','- **Persistent-line dominated coherence:** raw weighted coherence ≥ 0.3, at least 80% of cross-spectral weights in masked bins, and outside-line coherence <0.3.','', 'The catalog also retains zero-lag envelope r, signed peak r and lag, raw and masked spectral metrics, remaining SCM power fraction, peak coherent frequency, and longest contiguous run with coherence ≥0.6. Weights are sqrt(voltage PSD × SCM PSD), over 1–864 kHz. Masks are learned separately for each configuration.','','| Pair | Records | Median native r | Median envelope r | Median raw coherence | Median outside-line coherence |','|---|---:|---:|---:|---:|---:|']
 lines += [f"| {s['voltage']} / {s['scm']} | {s['n']} | {s['median_native_r0']:.3f} | {s['median_envelope_r0']:.3f} | {s['median_raw_coherence']:.3f} | {s['median_outside_line_coherence']:.3f} |" for s in summaries]
 lines+=['','Event-level flag counts (overlapping):','']+[f'- {k}: {v}' for k,v in diag['event_flag_counts'].items()]
 lines+=['','[Full pairwise coefficients](scm_voltage_correlations.csv) · [Thresholds, masks and diagnostics](scm_voltage_diagnostics.json)','']+[f'- {s}' for s in diag['caveats']]
 (args.out/'scm_voltage_report.md').write_text('\n'.join(lines));print(json.dumps({'event_flag_counts':diag['event_flag_counts'],'summaries':summaries}),flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser()
 for name in ['data','cache','out','figures']:ap.add_argument('--'+name,type=Path,required=True)
 run(ap.parse_args())
