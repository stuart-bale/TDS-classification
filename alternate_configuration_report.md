# Alternate V5 / SCM4 configuration

837 records, 819 unsaturated records fitted separately. Channels: V2, V5, V1V2, V3V4, SCM4.

Separate Ward fit with robust scaling and equal per-feature weights, k=2..4 with >=20 records per group. Forty record/feature resamples only if a partition qualifies. A numbering is independent of G.

Selected 2 groups; silhouette 0.524; median resampling ARI 0.988.

| A group | Records | Core | Medoid event |
|---|---:|---:|---:|
| A0 | 759 | 756 | 4363 |
| A1 | 60 | 60 | 4486 |

## Electron–voltage correlations

Within-burst analysis uses 511 unsaturated bursts with ≥32 counts. These coefficients use RMS envelopes in 128-sample bins and the sweep-phase baseline learned on other days with the same SWEAP status. No physical labels or calibration corrections are applied.

| Voltage | Pooled total / RMS Spearman | Median phase-residual envelope r at zero lag |
|---|---:|---:|
| V2 | -0.002 | -0.011 |
| V5 | 0.014 | 0.003 |
| V1V2 | 0.005 | -0.003 |
| V3V4 | -0.000 | -0.002 |

Pooled associations are sensitive to date, distance and instrument state. Daily coefficients and per-event lag correlations are in the diagnostic JSON and alternate_voltage_correlations.csv. These groups are exploratory: a small sample cannot establish stable physical populations. The primary G fit excludes this configuration; the separate A fit retains it without treating V5 as V4 or SCM4 as SCM5.