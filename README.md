# Reproduce the TDS atlas

Python 3.10+; dependencies: numpy, scipy, matplotlib, cdflib. The run used cdflib 1.3.12. Inputs remain read-only.

Run the scripts in this order (the cache is intermediate, not a deliverable):

```sh
python analyze_tds.py --data "/Users/bale/Dropbox/projects/PSP/Work/TDS/data" --cache ./cache --out .
python build_atlas.py --data "/Users/bale/Dropbox/projects/PSP/Work/TDS/data" --cache ./cache --out . --figures "/Users/bale/Dropbox/projects/PSP/Work/TDS/figures/waveform_classification"
python check_sensitivity.py --data "/Users/bale/Dropbox/projects/PSP/Work/TDS/data" --cache ./cache --out .
python add_context.py --data "/Users/bale/Dropbox/projects/PSP/Work/TDS/data" --cache ./cache --out . --figures "/Users/bale/Dropbox/projects/PSP/Work/TDS/figures/waveform_classification"
python alternate_configuration.py --data "/Users/bale/Dropbox/projects/PSP/Work/TDS/data" --cache ./cache --out . --figures "/Users/bale/Dropbox/projects/PSP/Work/TDS/figures/waveform_classification"
python electron_voltage.py --data "/Users/bale/Dropbox/projects/PSP/Work/TDS/data" --cache ./cache --out . --figures "/Users/bale/Dropbox/projects/PSP/Work/TDS/figures/waveform_classification"
python scm_voltage.py --data "/Users/bale/Dropbox/projects/PSP/Work/TDS/data" --cache ./cache --out . --figures "/Users/bale/Dropbox/projects/PSP/Work/TDS/figures/waveform_classification"
python publish_atlas.py --data "/Users/bale/Dropbox/projects/PSP/Work/TDS/data" --cache ./cache --out . --figures "/Users/bale/Dropbox/projects/PSP/Work/TDS/figures/waveform_classification"
```

Open atlas.html in a current Chrome, Edge, Safari, or Firefox with DecompressionStream support. The compressed embedded previews are unpacked locally; no network is used. If local storage is unavailable, labels remain in memory until you export them. Export annotations before closing. Catalog labels are provisional; see report.md for all thresholds, diagnostics, and limitations.

Group numbers are specific to each fitted dataset. Analog group numbers describe the current fit only; inspect medoids and feature medians. They carry no physical labels. Count group names remain generic pending expert review. Configuration exceptions retain actual channel labels and count analysis, but have no primary G or joint assignment; complete alternate records have a separate morphology analysis; see alternate_configuration_report.md; their primary canonical feature fields are NaN. Actual alternate-channel features and correlations are provided in separate tables. Missing channel preview slots are internal zeros and are rendered as absent, never as a measured zero signal. Do not carry physical interpretations from this run into a new dataset without review.

event_catalog.csv: one row per burst, full provenance and all group/status labels. features.csv: 89 analog features + 5 separate log-RMS amplitudes. count_features.csv: digital morphology and analog-count envelope relationships. phase_residual_features.csv: six features after phase conditioning. model_diagnostics.json: primary fit and resampling. counts_and_joint_diagnostics.json: count and joint partitions, metadata cross-tabs. sensitivity_diagnostics.json: learned spectral lines and cross-day phase profiles. source_inventory.json: input file sizes and SHA-256 hashes.
