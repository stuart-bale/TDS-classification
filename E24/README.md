# E24 ingestion and frozen-reference comparison

The input is the PI's completed June 2025 Dropbox delivery: 29 daily L2 TDS CDFs. June 4 has no file. Directory coverage does not certify completeness of spacecraft downlink.

All native analog samples were read to calculate the same 89 features used for the original E23 analysis. Each available digital series was checked against its CDF total. File hashes are retained, and size/mtime were checked against a pre-extraction snapshot. E24 IDs are allocated through `event_ids.json`; existing E23 IDs and the IDs allocated in the first E24 batch are preserved when later files are appended.

For compatible reference configurations, E24 features use the existing E23 medians, scales, family weights, clipping, and PCA basis. Reconstructed E23 coordinates are checked against the existing catalog before applying that basis to E24. The displayed group is the nearest original E23 centroid; it is a descriptive transfer, not an independent E24 cluster fit or a validated confidence classification. Saturated events are flagged and can be hidden. Alternate/incomplete configurations are shown in the waveform explorer but receive no PCA coordinates.

The PI's PC1 > 3.1 dust rule is transferred provisionally and explicitly marked. These axes contain the original spectral features, including power-supply contributions. They are not the pending harmonic-cleaned combined waveform/spectral PCA.

The viewer contains full-burst min/max envelopes, a native-sample zoom, logarithmically binned PSD previews, and digital count bins. Waveforms remain in engineering units. Full native series remain in the source CDFs; only previews and analysis tables are published. MAG projection and SPAN energy mapping have not been added to E24.

Reproduction: run `analyze_tds.py --stage extract` with `--data` pointing to the E24 folder, `--cache work/tds_e24`, and `--out outputs/tds_waveform_atlas/E24`; then run `work/build_e24.py` from the project directory with the existing E23 cache/model present. The latter expects `work/e24_template.html` and the input transfer snapshot. The extraction automatically verifies and reuses cached files.
