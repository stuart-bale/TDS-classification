"""MAG-supported voltage projections. No effective-length calibration assumed."""
from pathlib import Path
import json,sys,numpy as np
from scipy.io import readsav
from cdflib import CDF,cdfepoch
sys.path.insert(0,'outputs/tds_waveform_atlas')
from analyze_tds import P as PREFIX
P=Path('outputs/tds_github_pages');O=Path('outputs/tds_waveform_atlas/mag_projection');O.mkdir(exist_ok=True)
def basis(b):
 b=b/np.linalg.norm(b);a=np.array([-b[1],b[0],0.]);n=np.linalg.norm(a)
 a=a/n if n>1e-10 else np.array([1.,0.,0.]);return b,a,np.cross(b,a)
def samples(t,mt,my):
 j=np.searchsorted(mt,t);ix=np.array([j-1,j]);
 if j==0 or j==len(mt):return None
 if np.any(abs(mt[ix]-t)>.35) or not np.isfinite(my[ix]).all() or np.any(abs(my[ix])>1e5):return None
 return ix
def selfcheck():
 b,a,c=basis(np.array([1.,2.,3.]));q=np.array([b,a,c]);assert np.allclose(q@q.T,np.eye(3));x=np.random.default_rng(3).normal(size=(3,100));y=q@x;assert np.allclose(np.sum(x*x,0),np.sum(y*y,0));assert samples(1.,np.array([0.,2.]),np.ones((2,3))) is None
 # Constant transverse magnitude must have nonzero vector RMS, even though its standard deviation is zero.
 x=np.array([[1.,-1.],[0.,0.]]);assert np.sqrt(np.mean(np.sum(x*x,0)))==1
selfcheck();print('Projection and gap checks passed',flush=True)
s=readsav('/Users/bale/Dropbox/projects/PSP/Work/TDS/data/E23/MAG/MAG_SC_E23.tplot',python_dict=True)
for q in s['dq']:
 if q['name'].decode()=='psp_fld_l2_mag_SC_4_Sa_per_Cyc':mt=np.asarray(q['dh'][0]['x'],float);my=np.asarray(q['dh'][0]['y'],float).T;break
assert my.shape==(len(mt),3) and np.all(np.diff(mt)>0)
rows=json.loads((P/'pooled_waveform_points.json').read_text());out=[];last=None;c=None;g=np.deg2rad(40.);d=np.deg2rad(55.)
for i,r in enumerate(rows):
 base=dict(event=r['event'],encounter=r['encounter'])
 if r['encounter']!='E23':out.append(dict(**base,status='no_MAG_file'));continue
 if not all(ch in r['channels'] for ch in ['V1V2','V3V4']):out.append(dict(**base,status='missing_differential_channel'));continue
 if r['file']!=last:c=CDF(Path('/Users/bale/Dropbox/projects/PSP/Work/TDS/data/E23')/r['file']);last=r['file']
 ep=c.varget('epoch',startrec=r['record'],endrec=r['record'])[0];t=float(np.asarray(cdfepoch.unixtime([ep])).reshape(-1)[0])+r['n']/r['fs']/2;ix=samples(t,mt,my)
 if ix is None:out.append(dict(**base,status='no_close_bracketing_MAG'));continue
 B=my[ix].mean(0)
 if np.linalg.norm(B)<=0:out.append(dict(**base,status='zero_MAG'));continue
 bh,e1,e2=basis(B);x={};correction={}
 for ch in ['V1V2','V3V4','V2','V5']:
  if ch not in r['channels']:continue
  var=PREFIX+'Burst_Time_Series_'+ch+'_Engineering_mV';v=c.varget(var,startrec=r['record'],endrec=r['record']).ravel()[:r['n']].astype(float);assert len(v)==r['n'] and np.isfinite(v).all() and not np.any(v==c.varattsget(var).get('FILLVAL',-1e31));correction[ch]=float(np.median(v)-np.mean(v));x[ch]=v-np.mean(v)
 ex=x['V3V4']*np.cos(g)-x['V1V2']*np.cos(d);ey=-x['V3V4']*np.sin(g)-x['V1V2']*np.sin(d);full='V5' in x and 'V2' in x;ez=x['V5']-x['V2']-.5*x['V1V2'] if full else np.zeros(r['n']);v=np.array([ex,ey,ez]);par=bh@v;p1=e1@v;p2=e2@v
 rms=lambda z:float(np.sqrt(np.mean(z*z)))
 result=dict(**base,status='three_component_voltage_estimate' if full else 'xy_only_partial_parallel',B_SC_nT=B.tolist(),B_nT=float(np.linalg.norm(B)),Bhat=bh.tolist(),perp1_basis=e1.tolist(),perp2_basis=e2.tolist(),mag_offsets_s=(mt[ix]-t).tolist(),mag_direction_difference_deg=float(np.rad2deg(np.arccos(np.clip(np.dot(my[ix[0]],my[ix[1]])/np.linalg.norm(my[ix[0]])/np.linalg.norm(my[ix[1]]),-1,1)))),mean_correction_mV=correction,parallel_or_partial_rms_mV=rms(par),perp1_rms_mV=rms(p1),saturated=r['saturated'],fce_Hz=float(np.linalg.norm(B)*27.99249))
 if full:result.update(perp2_rms_mV=rms(p2),perp_vector_rms_mV=float(np.sqrt(np.mean(p1*p1+p2*p2))),total_rms_mV=float(np.sqrt(np.mean(np.sum(v*v,axis=0)))))
 if full:assert np.isclose(result['total_rms_mV']**2,rms(par)**2+result['perp_vector_rms_mV']**2)
 out.append(result)
 if i%500==0:print('MAG projections',i,flush=True)
from collections import Counter
summary=dict(total_events=len(rows),status_counts=dict(Counter(r['status'] for r in out)),method='Arithmetic average of the immediately bracketing MAG spacecraft-vector samples, each within 0.35 s of TDS epoch plus half burst duration. No distant nearest-sample fallback. Input voltages are mean-removed native engineering mV. PI rotation: X=V34 cos40-V12 cos55; Y=-V34 sin40-V12 sin55. When available, Z=V5-V2-V12/2. Parallel projection uses Bhat; one signed transverse axis is (-By,Bx,0)/sqrt(Bx^2+By^2), with X as fallback for axial B. Second transverse axis is Bhat cross first axis. Perpendicular-vector RMS is sqrt(mean(perp1^2+perp2^2)), only reported when Z is available. Without Z, parallel is explicitly partial and only the known in-plane transverse component is reported. These are voltage projections, not calibrated E in mV/m; differing antenna effective lengths and frequency responses are uncorrected. B averaging does not resolve magnetic variations within the 17 ms waveform. fce uses electron charge/mass and |B|; no density-based lower-hybrid estimate is supplied.')
data=dict(summary=summary,records=out);text=json.dumps(data,allow_nan=False,separators=(',',':'));(O/'context.json').write_text(text);(P/'mag_projection_context.json').write_text(text);print(json.dumps(summary,indent=2),flush=True)
