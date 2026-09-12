# Electron counts and voltage associations

All 16 files: 4,116 bursts. Across-burst comparisons use 3,845 unsaturated records with the common analog configuration. Within-burst comparisons require at least 32 counts and a matching-status phase baseline (2,992 bursts).

| Voltage | Pooled count-total / RMS Spearman | Adjusted rank correlation |
|---|---:|---:|
| V2 | 0.220 | 0.091 |
| V4 | 0.196 | 0.070 |
| V1V2 | 0.225 | 0.080 |
| V3V4 | 0.198 | 0.078 |

Rank residual correlation after day, SWEAP status, TWTA flag, selection type, radius rank and cubic spline of sweep phase (five interior knots). Descriptive observational adjustment, not a causal model.

| Voltage measure | Median raw zero-lag r | Median phase-residual zero-lag r | Number with residual absolute r ≥ 0.5 |
|---|---:|---:|---:|
| V2 envelope | 0.003 | 0.008 | 23 |
| V4 envelope | 0.001 | 0.005 | 18 |
| V1V2 envelope | -0.005 | 0.002 | 34 |
| V3V4 envelope | -0.006 | -0.000 | 33 |
| V2 signed_mean | -0.017 | -0.006 | 2 |
| V4 signed_mean | -0.012 | -0.005 | 4 |
| V1V2 signed_mean | -0.006 | -0.007 | 0 |
| V3V4 signed_mean | 0.001 | -0.000 | 3 |

Positive lag means counts follow voltage. The lag search covers ±1.067 ms in 66.7 µs steps. Voltage envelopes are RMS in 128-sample bins; signed means average those same samples. The high-frequency waveform phase is not resolved by this coarse signed trace.

Candidates below have ≥128 counts and are ranked by their largest absolute phase-residual correlation across four envelopes, four signed means and 33 lags. This selection inflates peak correlations; these are inspection candidates, not detections.

| Event | Counts | Channel / measure | Peak r | Counts lag (ms) | Matched donors |
|---|---:|---|---:|---:|---:|
| 127 | 17540 | V3V4 / envelope | -0.769 | 0.733 | 10 |
| 100 | 9456 | V1V2 / envelope | -0.752 | -0.133 | 11 |
| 2895 | 9553 | V2 / envelope | 0.738 | 0.867 | 8 |
| 2749 | 8600 | V1V2 / envelope | -0.735 | -1.000 | 10 |
| 2924 | 14329 | V1V2 / envelope | -0.726 | -1.000 | 9 |
| 993 | 8471 | V4 / envelope | 0.725 | 0.067 | 19 |
| 2927 | 14462 | V1V2 / envelope | -0.722 | 0.533 | 13 |
| 2866 | 13130 | V2 / envelope | 0.713 | -1.067 | 11 |

Matched donor comparisons were available for 10 events (at least 20 other bursts on the same day and SWEAP status, within 10 ms of sweep phase, with totals within a factor of four). The catalog retains zero-lag reference percentiles and exploratory circular-shift values. Neither establishes statistical significance under the instrument scan process.

The lag curves, distributions, daily comparisons and four candidate trace figures are saved in /Users/bale/Dropbox/projects/PSP/Work/TDS/figures/waveform_classification. The full per-burst coefficients are in electron_voltage_correlations.csv.

- SWEAP counts are scanned; phase residuals can retain scan structure.
- Correlation coefficients quantify association, not wave-particle coupling.
- Lag-peak selection searches eight traces and 33 lags; peak r is upward selected.
- Circular shifts assume stationarity and matched donor profiles are exchangeable only approximately; neither is a calibrated physical-detection test.
- Signed voltage uses 128-sample boxcar means, not resolved high-frequency voltage phase.
- No timing offset or transfer-function correction is applied; lag accuracy is limited by binning and instrument timing.

[Separate V5 / SCM4 analysis](alternate_configuration_report.md) covers the other analog configuration.

## Slow-trend sensitivity

At the originally selected trace and lag, removing a linear trend from both voltage and count residual gives:

| Event | Original peak r | Detrended r at same trace/lag | First-difference zero-lag r, same trace |
|---|---:|---:|---:|
| 127 | -0.769 | -0.510 | 0.071 |
| 100 | -0.752 | -0.023 | -0.008 |
| 2895 | 0.738 | 0.240 | 0.074 |
| 2749 | -0.735 | -0.326 | 0.057 |
| 2924 | -0.726 | -0.358 | 0.051 |
| 993 | 0.725 | 0.303 | -0.104 |
| 2927 | -0.722 | -0.559 | -0.053 |
| 2866 | 0.713 | 0.115 | 0.028 |

Differencing suppresses slow changes and amplifies counting noise; it probes timescale sensitivity, not whether the original association is physical.

## Native-cadence check

The catalog also contains direct signed voltage/count correlations at 1.92 MSa/s and peak magnitude-squared Welch coherence between 1 and 864 kHz. These use original engineering voltage samples and digital count increments. Spectra use 4096-sample windows with 50% overlap; frequency-response phase corrections and instrumental relative timing are unavailable. Peak coherence searches many bins and is not a detection statistic.

| Voltage | Median native zero-lag r | Largest absolute native r | Median peak coherence |
|---|---:|---:|---:|
| V2 | -0.00107 | 0.134 | 0.444 |
| V4 | -0.00044 | 0.146 | 0.446 |
| V1V2 | -0.00013 | 0.139 | 0.446 |
| V3V4 | 0.00014 | 0.123 | 0.446 |