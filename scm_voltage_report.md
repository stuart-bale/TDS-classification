# SCM–voltage association flags

4,664 of 4,715 bursts assessed (18,656 voltage/SCM pairs). Saturation and missing-SCM records remain unassessed. All four available voltages are tested against SCM5 or SCM4, keeping configurations separate.

Flags are deliberately descriptive and may overlap:

- **Waveform correlation:** absolute full-sample zero-lag Pearson r ≥ 0.5.
- **Envelope association:** absolute peak RMS-envelope r ≥ 0.5 across ±1.067 ms (66.7 µs bins); positive lag means SCM follows voltage.
- **Coherence outside persistent lines:** PSD-weighted mean magnitude-squared coherence ≥ 0.3 after excluding learned narrow-line bins. This does not necessarily mean broadband coherence.
- **Persistent-line dominated coherence:** raw weighted coherence ≥ 0.3, at least 80% of cross-spectral weights in masked bins, and outside-line coherence <0.3.

The catalog also retains zero-lag envelope r, signed peak r and lag, raw and masked spectral metrics, remaining SCM power fraction, peak coherent frequency, and longest contiguous run with coherence ≥0.6. Weights are sqrt(voltage PSD × SCM PSD), over 1–864 kHz. Masks are learned separately for each configuration.

| Pair | Records | Median native r | Median envelope r | Median raw coherence | Median outside-line coherence |
|---|---:|---:|---:|---:|---:|
| V2 / SCM5 | 3845 | -0.092 | -0.000 | 0.815 | 0.072 |
| V4 / SCM5 | 3845 | -0.021 | 0.002 | 0.510 | 0.072 |
| V1V2 / SCM5 | 3845 | 0.020 | -0.002 | 0.608 | 0.072 |
| V3V4 / SCM5 | 3845 | -0.001 | 0.001 | 0.663 | 0.073 |
| V2 / SCM4 | 819 | -0.208 | -0.001 | 0.882 | 0.071 |
| V5 / SCM4 | 819 | -0.095 | -0.004 | 0.934 | 0.071 |
| V1V2 / SCM4 | 819 | 0.092 | -0.002 | 0.783 | 0.071 |
| V3V4 / SCM4 | 819 | -0.006 | -0.007 | 0.775 | 0.069 |

Event-level flag counts (overlapping):

- waveform correlation: 174
- envelope association: 43
- coherence outside persistent lines: 174
- persistent-line dominated coherence: 3010
- no threshold crossed: 1462
- not assessed: 51

[Full pairwise coefficients](scm_voltage_correlations.csv) · [Thresholds, masks and diagnostics](scm_voltage_diagnostics.json)

- Flags are descriptive screening thresholds, not significance tests or physical wave identifications.
- SCM4 and SCM5 are evaluated separately; saturated or missing-SCM records are not assessed.
- Persistent lines are learned separately per configuration (>8 dB local excess in >=60% of unsaturated records, +/-1 Welch bin); persistence does not establish instrumental origin.
- Welch coherence uses 15 overlapping 4096-sample segments; individual spectral maxima are upward selected and are not independently significant.
- Magnitude-squared coherence is averaged with sqrt(voltage PSD × SCM PSD) weights; removing lines changes both bins and weights.
- No transfer-function phase correction or instrument-relative timing correction is applied. Coherence does not establish physical E/B polarization, propagation, impedance or Poynting flux.