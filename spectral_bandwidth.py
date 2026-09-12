"""Per-channel Welch peak and local-background FWHM; no PCA refit."""
import json,csv
from pathlib import Path
import numpy as np
from scipy import signal
from cdflib import CDF
from analyze_tds import P,channel_names,waveform
OUT=Path(__file__).resolve().parent

def measure(x,fs,limit):
 f,s=signal.welch(x,fs,window='hann',nperseg=4096,noverlap=2048,detrend='constant')
 keep=(f>0)&(f<=limit);f=f[keep];s=s[keep];df=fs/4096
 if not np.any(s>0):return {'flags':['zero_power']}
 k=int(np.argmax(s));radius=max(20,int(np.ceil(.2*f[k]/df)));idx=np.arange(max(0,k-radius),min(len(s),k+radius+1));bgidx=idx[abs(idx-k)>3];bg=float(np.median(s[bgidx]));level=bg+(s[k]-bg)/2
 l=k;r=k
 while l>0 and s[l]>level:l-=1
 while r<len(s)-1 and s[r]>level:r+=1
 flags=[]
 if s[l]>level or s[r]>level:width=None;left=None;right=None;fraction=None;flags.append('band_edge_censored')
 else:
  left=float(f[l]+(level-s[l])/(s[l+1]-s[l])*df);right=float(f[r-1]+(level-s[r-1])/(s[r]-s[r-1])*df);width=right-left
  ff=np.r_[left,f[(f>left)&(f<right)],right];ss=np.interp(ff,f,s);fraction=float(np.trapz(ss,ff)/(s.sum()*df))
  if width<=2*df:flags.append('resolution_limited')
 peaks,props=signal.find_peaks(s,prominence=max((s[k]-bg)*.25,0),distance=3)
 if len(peaks)>1:flags.append('multiple_prominent_peaks')
 if any(abs(f[k]-hz)<=df for hz in np.arange(150000,fs/2,50000)):flags.append('near_known_line_grid')
 return {'peak_hz':float(f[k]),'fwhm_hz':width,'fractional_bandwidth':None if width is None else width/f[k],'left_hz':left,'right_hz':right,'peak_power_fraction':fraction,'background_psd':bg,'peak_to_background':float(s[k]/max(bg,1e-300)),'frequency_bin_hz':df,'prominent_peak_count':len(peaks),'flags':flags}

def main():
 rows=json.loads(Path('work/tds/final_records.json').read_text());results=[];c=None;last=None
 for i,row in enumerate(rows):
  if row['file']!=last:c=CDF(Path('/Users/bale/Dropbox/projects/PSP/Work/TDS/data/E23')/row['file']);last=row['file']
  n=int(row['N_Samples_per_Channel']);rec=int(row['record']);names=channel_names(c,rec);x=waveform(c,rec,n)
  for j,ch in enumerate(names):
   if not ch:continue
   m=measure(x[j],row['Sample_Speed'],min(row['Low_Pass_Filter'],row['Sample_Speed']/2));m.update(event=row['event'],channel=ch,utc=row['utc'],file=row['file'],record=rec,pc1=row.get('pc1'),dust=row.get('group',-1)>=0 and row.get('pc1',-1)>3.1,saturated=bool(row['Burst_Saturation_Flag']))
   results.append(m)
  if i%500==0:print(i,'of',len(rows),flush=True)
 (OUT/'spectral_bandwidth.json').write_text(json.dumps(results,allow_nan=False,separators=(',',':')))
 keys=list(dict.fromkeys(k for r in results for k in r))
 with (OUT/'spectral_bandwidth.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({**r,'flags':'; '.join(r['flags'])} for r in results)
 print('DONE',len(results),flush=True)
if __name__=='__main__':main()
