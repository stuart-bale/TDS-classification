# Alternate V5 / SCM4 configuration

241 records, 234 unsaturated records fitted separately. Channels: V2, V5, V1V2, V3V4, SCM4.

Separate Ward fit with robust scaling and equal per-feature weights, k=2..4 with >=20 records per group. Forty record/feature resamples only if a partition qualifies. A numbering is independent of G.

No candidate partition met the minimum 20 records per group. A0 is an unpartitioned browsing pool, not a discovered class.

| A group | Records | Core | Medoid event |
|---|---:|---:|---:|
| A0 | 234 | 0 | 3099 |

## Electron–voltage correlations

Within-burst analysis uses 151 unsaturated bursts with ≥32 counts. These coefficients use RMS envelopes in 128-sample bins and the sweep-phase baseline learned on other days with the same SWEAP status. No physical labels or calibration corrections are applied.

| Voltage | Pooled total / RMS Spearman | Median phase-residual envelope r at zero lag |
|---|---:|---:|
| V2 | 0.119 | -0.013 |
| V5 | 0.047 | 0.006 |
| V1V2 | 0.104 | -0.013 |
| V3V4 | 0.054 | -0.009 |

Pooled associations are sensitive to date, distance and instrument state. Daily coefficients and per-event lag correlations are in the diagnostic JSON and alternate_voltage_correlations.csv. These groups are exploratory: a small sample cannot establish stable physical populations. The primary G fit excludes this configuration; the separate A fit retains it without treating V5 as V4 or SCM4 as SCM5.