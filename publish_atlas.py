"""Assemble a standalone local HTML explorer, report, and reproducibility notes."""
import argparse,base64,gzip,json,html,sys,csv
from pathlib import Path
import numpy as np
from collections import Counter
from analyze_tds import dump,transform
from build_atlas import plot_event,savecsv
from build_methods import build_methods
from encounter_tags import tag_encounter
from expert_labels import apply_expert_labels

NAMES=['Sustained oscillatory / spectrally concentrated',
       'Irregular electric fluctuations / broader spectra',
       'Mixed oscillations / broader SCM spectrum']
CNAMES=['Counts spread through the record','Counts concentrated in part of the record']

def publish(args):
    rows=json.loads((args.cache/'final_records.json').read_text());prev=json.loads((args.cache/'previews.json').read_text())
    model=json.loads((args.out/'model_diagnostics.json').read_text());cj=json.loads((args.out/'counts_and_joint_diagnostics.json').read_text())
    sens=json.loads((args.out/'sensitivity_diagnostics.json').read_text());z=np.load(args.cache/'features.npz');keys=z['keys'].tolist()
    alternate=json.loads((args.out/'alternate_configuration_diagnostics.json').read_text())
    altnote=(f"The separate V5 / SCM4 fit selected {alternate['k']} provisional A groups." if alternate['accepted_partition'] else 'No candidate met the minimum group-size rule; A0 is an unpartitioned browsing pool.')
    global NAMES,CNAMES
    NAMES=[f'Morphology group {i}' for i in range(model['k'])]
    CNAMES=[f'Count-shape group {i}' for i in range(cj['counts']['k'])]
    impulse=10**np.mean(z['X'][:,[keys.index(ch+'__log_crest') for ch in ['V2','V4','V1V2','V3V4']]],axis=1)
    for i,r in enumerate(rows):
        r['impulsiveness']=float(impulse[i]) if np.isfinite(impulse[i]) else None;r['group_description']=NAMES[r['group']] if r['group']>=0 else altnote if r.get('alternate_group',-1)>=0 else 'Unassigned: flagged or incomplete analog configuration'
        r['count_group_description']=CNAMES[r['count_group']] if r['count_group']>=0 else r['count_status']
    alternate=json.loads((args.out/'alternate_configuration_diagnostics.json').read_text())
    correlations=json.loads((args.out/'electron_voltage_diagnostics.json').read_text())
    for r in rows:r['voltage_count_peak_r']=None;r['voltage_count_peak_abs_r']=None;r['voltage_count_peak_lag_ms']=None;r['voltage_count_peak_trace']='Unavailable'
    for v in csv.DictReader((args.out/'electron_voltage_correlations.csv').open()):
        rows[int(v['event'])].update(voltage_count_peak_r=float(v['residual_r_at_peak']),voltage_count_peak_abs_r=float(v['max_residual_abs_r']),voltage_count_peak_lag_ms=float(v['counts_lag_ms']),voltage_count_peak_trace=v['peak_channel']+' '+v['peak_measure'])
    scm=json.loads((args.out/'scm_voltage_diagnostics.json').read_text())
    for row in rows:
        tag_encounter(row)
        row.setdefault('counts_present',True)
        row.setdefault('Burst_Time_Series_SWEAP_Exists_Flag',1)
    apply_expert_labels(rows)
    savecsv(args.out/'event_catalog.csv',rows)
    selected=set(json.loads((args.cache/'selected_examples.json').read_text()))
    selected.update(g['medoid'] for g in sens['phase_conditioned_counts']['groups'])
    selected.update(int(i) for i in np.argsort(np.nan_to_num(impulse,nan=-1))[-12:])
    for old in (args.figures/'examples').glob('event_*.png'):
        if int(old.stem.split('_')[-1]) not in selected:old.unlink()
    if not args.skip_figures:
        for i in sorted(selected):plot_event(args,rows[i],args.figures/'examples'/f'event_{i:04d}.png')
    examples={i:(args.figures/'examples'/f'event_{i:04d}.png').as_uri() for i in selected}
    packed=[]
    for i,p in enumerate(prev):
        lo=np.array(p['lo']).reshape(5,256,2).min(axis=2)
        hi=np.array(p['hi']).reshape(5,256,2).max(axis=2)
        cp=p['counts'].copy();cp['zoom']=cp['zoom'][256:768]
        packed.append(dict(channels=p.get('channels',['V2','V4','V1V2','V3V4','SCM5']),lo=lo.tolist(),hi=hi.tolist(),zoom=np.array(p['zoom'])[:,256:768].tolist(),
             zoom_start=p['zoom_start']+256,freq=p['freq'],psd_db=p['psd_db'],rms=p['rms'],counts=cp))
    weighted,_=transform(z['X'],np.array(model['scaling_median']),np.array(model['scaling_iqr']),keys)
    scm_pairs={}
    for v in csv.DictReader((args.out/'scm_voltage_correlations.csv').open()):
        pair={k:v[k] for k in ['voltage','scm','flags']}
        pair.update({k:float(v[k]) for k in ['native_r0','envelope_r0','envelope_peak_r','scm_lag_ms','weighted_coherence_raw','weighted_coherence_outside_lines','scm_power_remaining_fraction','peak_coherence_frequency_hz']})
        scm_pairs.setdefault(int(v['event']),[]).append(pair)
    data={'scm_pairs':scm_pairs,'rows':rows,'p':packed,'features':np.round(weighted,7).tolist(),'model':model,'counts':cj,'sensitivity':sens,'names':NAMES,'cnames':CNAMES,'examples':examples,'alternate':alternate,'scm':scm}
    encoded=base64.b64encode(gzip.compress(json.dumps(data,separators=(',',':'),allow_nan=False).encode(),compresslevel=6)).decode()
    overview=base64.b64encode((args.figures/'overview.png').read_bytes()).decode()
    template=(Path(__file__).parent/'atlas_template.html.in').read_text()
    radial=base64.b64encode((args.figures/'radial_context.png').read_bytes()).decode()
    scmfig=base64.b64encode((args.figures/'scm_voltage_associations.png').read_bytes()).decode()
    scmcards=' '.join(f'<button data-pick="{r["event"]}">Event {r["event"]} · {r["voltage"]}/{r["scm"]}</button>' for r in scm['candidate_events'])
    scmdetail=''
    for r in scm['candidate_events']:
        event=r['event'];png=args.figures/f'scm_voltage_event_{event:04d}.png';encoded_png=base64.b64encode(png.read_bytes()).decode()
        scmdetail+=f'<details><summary>SCM association: event {event}</summary><img class="overview" src="data:image/png;base64,{encoded_png}" alt="SCM and voltage waveform, envelope and coherence comparison"></details>'
    template=template.replace('__SCMFIG__',scmfig).replace('__SCMCARDS__',scmcards).replace('__SCMDETAIL__',scmdetail)
    corrfig=base64.b64encode((args.figures/'electron_voltage_correlations.png').read_bytes()).decode()
    cards=''.join(f'<button data-pick="{r["event"]}">Event {r["event"]} · r = {r["residual_r_at_peak"]:.2f}</button> ' for r in correlations['candidate_events'])
    detail=''
    for r in correlations['candidate_events'][:4]:
        event=r['event'];png=args.figures/f'electron_voltage_event_{event:04d}.png'
        encoded_png=base64.b64encode(png.read_bytes()).decode()
        detail+=f'<details><summary>Correlation traces: event {event}</summary><img class="overview" src="data:image/png;base64,{encoded_png}" alt="Count and voltage-envelope traces"></details>'

    altfig=base64.b64encode((args.figures/'alternate_voltage_correlations.png').read_bytes()).decode()
    altcards=' '.join(f'<button data-pick="{g["medoid"]}">A{g["group"]} · {g["n"]} bursts · medoid {g["medoid"]}</button>' for g in alternate['groups'])
    ar=correlations['across_burst'];pool=[v['pooled_spearman'] for v in ar];adjusted=[v['adjusted_rank_correlation'] for v in ar]
    template=template.replace('__POOLR__',f'{min(pool):.2f}–{max(pool):.2f}').replace('__ADJR__',f'{min(adjusted):.2f}–{max(adjusted):.2f}')
    template=template.replace('__ALTFIG__',altfig).replace('__ALTCARDS__',altcards)
    template=template.replace('__CORRFIG__',corrfig).replace('__CORRCARDS__',cards).replace('__CORRDETAIL__',detail)
    template=template.replace('__PAYLOAD__',encoded).replace('__OVERVIEW__',overview).replace('__RADIALFIG__',radial)
    st=Counter(r['status'] for r in rows);cd=cj['counts'];rd=sens['phase_conditioned_counts']
    group_table='\n'.join(f"| G{g['group']} | {NAMES[g['group']]} | {g['n']} | {g['core']} | {g['medoid']} |" for g in model['groups'])
    alt='\n'.join(f"| {k.replace('_',' ')} | {v['ARI']:.3f} |" for k,v in model['sensitivity'].items())
    counts_table='\n'.join(f"| C{g['group']} | {CNAMES[g['group']]} | {g['n']} | {g['core']} | {g['medoid']} |" for g in cd['groups'])
    residual_table='\n'.join(f"| R{g['group']} | {g['n']} | {g['core']} | {g['median_residual_rms']:.3f} | {g['medoid']} |" for g in rd['groups'])
    refitnote='No prior-run comparison is available for this extraction.'
    if (args.out/'refit_comparison.json').exists():
        comparison=json.loads((args.out/'refit_comparison.json').read_text())
        if comparison['current_records']==len(rows):
            refitnote=f"Relative to the previous {comparison['previous_records']:,}-record run, the selected group count changes from {comparison['old_k']} to {comparison['new_k']}. Agreement on the same {comparison['n_compared']:,} previously fitted records is ARI {comparison['ARI_on_same_preexisting_fit_records']:.3f}. Group numbering and membership must be reviewed after refitting; event IDs remain stable."
    replacements={'__REFITNOTE__':refitnote,'__ALTNOTE__':altnote,
       '__N__':str(len(rows)), '__DATES__':min(r['utc'] for r in rows)[:10]+' to '+max(r['utc'] for r in rows)[:10],
       '__K__':str(model['k']), '__CORE__':str(st['core']), '__AMBIG__':str(st['ambiguous']), '__OUTLIERS__':str(st['outlier']), '__SAT__':str(st['saturated']),
       '__SIL__':f"{model['silhouette']:.3f}", '__ARI__':f"{model['stability_ARI']['median']:.3f}",
       '__LINEPOWER__':f"{sens['persistent_lines']['channels'][4]['median_removed_power_fraction']:.1%}",
       '__LINEARI__':f"{sens['persistent_lines']['fixed_k_ARI_vs_raw']:.3f}",
       '__R2__':f"{rd['cross_day_log1p_total_R2']:.3f}", '__GROUPTABLE__':group_table, '__COUNTTABLE__':counts_table,
       '__RESTABLE__':residual_table,'__ALTTABLE__':alt, '__ZEROS__':str(cd['status_counts']['zero counts']), '__SPARSE__':str(cd['status_counts']['too few counts']),
       '__COUNTN__':str(sum(g['n'] for g in cd['groups'])), '__COUNTSIL__':f"{cd['silhouette']:.3f}", '__COUNTARI__':f"{cd['resample_ARI_median']:.3f}",
       '__RATEARI__':f"{cd['ARI_with_count_rate']:.3f}", '__RESARI__':f"{rd['resample_ARI_median']:.3f}",
       '__JOINTN__':str(cj['joint']['n']), '__JOINTARI__':f"{cj['joint']['ARI_vs_analog']:.3f}", '__CROSSARI__':f"{cj['joint']['counts_vs_analog_ARI']:.3f}",
       '__FIGURES__':str(args.figures), '__DAYCOUNTS__':', '.join(f"{day}: {n}" for day,n in sorted(Counter(r['day'] for r in rows).items())),
       '__BESTK__':str(max(model['candidates'],key=lambda v:v['silhouette'])['k']),
       '__AVGCOUNTS__':str(model['sensitivity']['average_linkage']['sizes'])}
    replacements.update(__RANGE__=f"{min(r['solar_distance_Rs'] for r in rows):.2f}–{max(r['solar_distance_Rs'] for r in rows):.2f}",__TWTA__=str(sum(bool(r['SC_TWTA_On_Flag']) for r in rows)),__SATFLAGS__=str(sum(bool(r['Burst_Saturation_Flag']) for r in rows)),__CONFIG__=str(sum(r['group']<0 for r in rows)),__C0N__=str(cd['groups'][0]['n']),__C1N__=str(cd['groups'][1]['n']),
        __R0N__=str(rd['groups'][0]['n']),__R1N__=str(rd['groups'][1]['n']),
        __PCAPERCENT__=f"{sum(model['pca_explained']):.0%}",
        __DAYOPTIONS__=''.join(f'<option value="{d}">{d[6:]} March</option>' for d in sorted({r['day'] for r in rows})))
    for key,value in replacements.items():template=template.replace(key,value)
    (args.out/'atlas.html').write_text(template)
    radial_info=json.loads((args.out/'radial_context.json').read_text())
    replacements['__RADIUSTABLE__']='\n'.join(f"| {f['day']} | {f['min_Rs']:.2f}–{f['max_Rs']:.2f} | {f['directions']} |" for f in radial_info['files'])
    report=(Path(__file__).parent/'report_template.md').read_text()
    for key,value in replacements.items():report=report.replace(key,value)
    (args.out/'report.md').write_text(report)
    build_methods(args.out)
    readme=f'''# Reproduce the TDS atlas

Python 3.10+; dependencies: numpy, scipy, matplotlib, cdflib. The run used cdflib 1.3.12. Inputs remain read-only.

Run the scripts in this order (the cache is intermediate, not a deliverable):

```sh
python analyze_tds.py --data "{args.data}" --cache ./cache --out .
python build_atlas.py --data "{args.data}" --cache ./cache --out . --figures "{args.figures}"
python check_sensitivity.py --data "{args.data}" --cache ./cache --out .
python add_context.py --data "{args.data}" --cache ./cache --out . --figures "{args.figures}"
python alternate_configuration.py --data "{args.data}" --cache ./cache --out . --figures "{args.figures}"
python electron_voltage.py --data "{args.data}" --cache ./cache --out . --figures "{args.figures}"
python scm_voltage.py --data "{args.data}" --cache ./cache --out . --figures "{args.figures}"
python publish_atlas.py --data "{args.data}" --cache ./cache --out . --figures "{args.figures}"
```

Open atlas.html in a current Chrome, Edge, Safari, or Firefox with DecompressionStream support. The compressed embedded previews are unpacked locally; no network is used. If local storage is unavailable, labels remain in memory until you export them. Export annotations before closing. Catalog labels are provisional; see report.md for all thresholds, diagnostics, and limitations.

Group numbers are specific to each fitted dataset. Analog group numbers describe the current fit only; inspect medoids and feature medians. They carry no physical labels. Count group names remain generic pending expert review. Configuration exceptions retain actual channel labels and count analysis, but have no primary G or joint assignment; complete alternate records have a separate morphology analysis; see alternate_configuration_report.md; their primary canonical feature fields are NaN. Actual alternate-channel features and correlations are provided in separate tables. Missing channel preview slots are internal zeros and are rendered as absent, never as a measured zero signal. Do not carry physical interpretations from this run into a new dataset without review.

event_catalog.csv: one row per burst, full provenance and all group/status labels. features.csv: 89 analog features + 5 separate log-RMS amplitudes. count_features.csv: digital morphology and analog-count envelope relationships. phase_residual_features.csv: six features after phase conditioning. model_diagnostics.json: primary fit and resampling. counts_and_joint_diagnostics.json: count and joint partitions, metadata cross-tabs. sensitivity_diagnostics.json: learned spectral lines and cross-day phase profiles. source_inventory.json: input file sizes and SHA-256 hashes.
'''
    (args.out/'README.md').write_text(readme)
    import scipy,cdflib,matplotlib,platform
    dump(args.out/'runtime_versions.json',dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,cdflib=cdflib.__version__,matplotlib=matplotlib.__version__))
    print('Published atlas bytes',(args.out/'atlas.html').stat().st_size,'figure count',len(selected)+1,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);ap.add_argument('--cache',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--figures',type=Path,required=True);ap.add_argument('--skip-figures',action='store_true');publish(ap.parse_args())
