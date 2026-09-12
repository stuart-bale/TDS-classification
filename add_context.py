"""Attach CDF spacecraft-reported solar distance as withheld physical context."""
import argparse,json
from collections import Counter
from pathlib import Path
import numpy as np
from cdflib import CDF
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from analyze_tds import P,dump
from build_atlas import savecsv,COLORS

def run(args):
    rows=json.loads((args.cache/'final_records.json').read_text());summary=[]
    for file in sorted({r['file'] for r in rows}):
        c=CDF(args.data/file)
        km=c.varget(P+'SC_Solar_Distance');known=c.varget(P+'SC_Solar_Distance_Known_Flag');valid=c.varget(P+'SC_Solar_Distance_Valid_Flag');direction=c.varget(P+'SC_Solar_Distance_Direction')
        for r in rows:
            if r['file']!=file:continue
            i=r['record'];r['solar_distance_known']=int(known[i]);r['solar_distance_valid']=int(valid[i])
            r['solar_distance_km']=float(km[i]) if known[i]==1 and valid[i]==1 else None
            r['solar_distance_Rs']=float(km[i]/695700) if r['solar_distance_km'] is not None else None
            r['solar_distance_AU']=float(km[i]/149597870.7) if r['solar_distance_km'] is not None else None
            r['radial_direction']='inbound' if direction[i]==0 else 'outbound' if direction[i]==1 else 'unknown'
        rr=[r for r in rows if r['file']==file];v=[r['solar_distance_Rs'] for r in rr if r['solar_distance_Rs'] is not None]
        summary.append({'file':file,'day':rr[0]['day'],'n':len(rr),'min_Rs':min(v),'max_Rs':max(v),'directions':dict(Counter(r['radial_direction'] for r in rr))})
    dump(args.cache/'final_records.json',rows);savecsv(args.out/'event_catalog.csv',rows)
    valid=[r for r in rows if r['solar_distance_Rs'] is not None];closest=min(valid,key=lambda r:r['solar_distance_Rs'])
    diag={'source_variable':P+'SC_Solar_Distance','source_description':'Spacecraft-reported Sun distance at each burst, from CDF telemetry. No SPICE state vectors loaded.',
          'solar_radius_km':695700,'AU_km':149597870.7,'valid_records':len(valid),'files':summary,
          'closest_record':{k:closest[k] for k in ['event','utc','solar_distance_km','solar_distance_Rs','radial_direction']},
          'used_for_clustering':False}
    dump(args.out/'radial_context.json',diag)
    fig,axs=plt.subplots(2,2,figsize=(13,8),layout='constrained');k=max(r['group'] for r in rows)+1
    for g in range(k):
        for direction,marker in [('inbound','v'),('outbound','o')]:
            rr=[r for r in valid if r['group']==g and r['radial_direction']==direction and r['status'] in ['core','ambiguous','outlier']]
            axs[0,0].scatter([np.datetime64(r['utc']) for r in rr],[r['solar_distance_Rs'] for r in rr],s=10,color=COLORS[g],marker=marker,alpha=.6,label=f'G{g} {direction}')
            axs[1,0].scatter([r['solar_distance_Rs'] for r in rr],[r['V3V4__rms'] for r in rr],s=10,color=COLORS[g],marker=marker,alpha=.5)
    excluded=[r for r in valid if r['group']<0]
    axs[0,0].scatter([np.datetime64(r['utc']) for r in excluded],[r['solar_distance_Rs'] for r in excluded],s=12,color='#888',marker='x',label='Alternate / incomplete')
    axs[0,0].xaxis.set_major_formatter(mdates.DateFormatter('%d Mar\n%H:%M'));axs[0,0].set(ylabel='Distance from Sun center (R☉)',title='Distance at each recorded burst');axs[0,0].legend(fontsize=7,ncol=2)
    axs[1,0].set(xlabel='Heliocentric distance (R☉)',ylabel='V3−V4 engineering RMS (mV)',yscale='log',title='Amplitude and morphology versus distance')
    edges=np.array([9.5,12,15,20,30,40,60,80,100,120]);data=[r for r in valid if r['status'] in ['core','ambiguous','outlier']];tab=np.array([[sum(r['group']==g and a<=r['solar_distance_Rs']<b for r in data) for a,b in zip(edges[:-1],edges[1:])] for g in range(k)])
    totals=tab.sum(axis=0);bottom=np.zeros(len(totals))
    for g in range(k):
        height=tab[g]/np.maximum(totals,1);axs[0,1].bar(np.arange(len(totals)),height,bottom=bottom,color=COLORS[g],label=f'G{g}');bottom+=height
    axs[0,1].set_xticks(np.arange(len(totals)),[f'{a:g}–{b:g}\nn={n}' for a,b,n in zip(edges[:-1],edges[1:],totals)]);axs[0,1].tick_params(axis='x',labelsize=8);axs[0,1].set(xlabel='Distance bin (R☉)',ylabel='Fraction of recorded unsaturated bursts',title='Group mix · common configuration only');axs[0,1].legend(fontsize=8)
    for g in sorted(set(r['count_group'] for r in rows)):
        rr=[r for r in valid if r['count_group']==g];axs[1,1].scatter([r['solar_distance_Rs'] for r in rr],[r['counts_total']+1 for r in rr],s=10,color=COLORS[g] if g>=0 else '#999',alpha=.5,label=f'C{g}' if g>=0 else 'Too few / zero counts')
    axs[1,1].set(xlabel='Heliocentric distance (R☉)',ylabel='SPAN-e total counts + 1',yscale='log',title='Digital counts: distance and scan phase both matter');axs[1,1].legend(fontsize=8)
    for ax in axs.ravel():ax.grid(alpha=.15)
    fig.suptitle(f"PSP / TDS · Radial context for {len(rows):,} bursts\nSpacecraft-reported CDF distance; R☉ = 695,700 km. Distance excluded from clustering.",fontsize=13)
    fig.savefig(args.figures/'radial_context.png',dpi=160);plt.close(fig)
    print(json.dumps(diag),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);ap.add_argument('--cache',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--figures',type=Path,required=True);run(ap.parse_args())
