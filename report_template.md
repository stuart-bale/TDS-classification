# TDS waveform classification: exploratory atlas

**__N__ six-channel bursts, __DATES__.** Per-file record counts: __DAYCOUNTS__. This is a reproducible morphology atlas and event catalog, with analog, digital-count, and joint partitions. No plasma-wave or dust examples, templates, or physical labels were supplied to the algorithms. Group numbers describe the current fit only. Inspect medoids and feature medians; no physical origins are assigned.

## Findings

The expanded dataset favors **__K__ broad analog groups**. Mean silhouette is **__SIL__**, and median resampling adjusted Rand index (ARI) is **__ARI__**. There is useful recurring structure, but substantial overlap and representation dependence. This is **not a validated physical taxonomy**. The catalog retains **__CORE__ core, __AMBIG__ ambiguous, __OUTLIERS__ outlier, __SAT__ saturation-status, and __CONFIG__ configuration-exception** records. In total __SATFLAGS__ records have the saturation flag: some also have an alternate or incomplete configuration. These statuses are mutually exclusive; the original saturation flag is preserved.

The digital counts show a repeatable dependence on SWEAP sweep phase. A profile trained on the other days within each SWEAP status predicts log(1 + total counts) with pooled held-out-day **R² = __R2__**. This is a prediction score on log-transformed totals, not a fraction of individual electrons or raw count variance. Repeated rises and falls across sweep phase are visible in the overview; a monotonic correlation alone misses this structure.

Detected persistent spectral-line bins contain a median **__LINEPOWER__ of SCM5 Welch power**. The corresponding electric-channel fractions are much smaller. Masking these bins in spectral descriptors changes the analog partition substantially (ARI **__LINEARI__**). Persistence alone does not prove an instrumental origin. Original engineering waveforms are retained, and masking is only a sensitivity experiment.

The initial two-day fit favored five smaller groups. The three-day fit favored three broader groups; the expanded result above supersedes it. That change, and the weak independent-day agreement for some days, show that these boundaries are provisional and sample-dependent.

__REFITNOTE__

See refit_comparison.json for the old-to-new assignment table.

## Heliocentric distance

Every record has both solar-distance-known and solar-distance-valid flags set. The values below come directly from `PSP_FLD_L2_TDS_WF_SC_Solar_Distance`, described by the CDF as spacecraft-reported Sun distance at the burst time. They are not independently recomputed from SPICE. Conversions use R☉ = 695,700 km and AU = 149,597,870.7 km.

| Date | Distance from Sun center (R☉) | Direction counts |
|---|---:|---|
__RADIUSTABLE__

Distance spans approximately __RANGE__ R☉. The March 22 file includes the inward-to-outward transition; subsequent records are outbound. The atlas shows radius and direction for each burst, filters by both, and includes radial group/amplitude comparisons. These context variables were withheld from clustering. Differences between days can therefore be examined as possible radial trends, while accounting for selection, instrument effects, and overlapping inbound/outbound coverage. The closest recorded burst is not necessarily the exact perihelion instant.

Full position and velocity vectors have not been loaded; radial context here uses the per-burst CDF distance.

## Alternate analog configuration

[Separate V5 / SCM4 results](alternate_configuration_report.md) include a separate morphology assessment and voltage/count correlations. __ALTNOTE__ The primary fit below retains a consistent channel configuration.

## Analog groups

Parent-group counts exclude saturated and configuration-exception records but include ambiguous and outlier members. “Core” means stable enough for the stated exploratory policy, not physically confirmed. Saturated records have a nearest-centroid suggestion only.

| Group | Descriptive morphology | Parent-group records | Core records | Medoid event |
|---|---|---:|---:|---:|
__GROUPTABLE__

Group numbering belongs to this expanded fit and is independent of earlier group numbering. Inspect the listed medoids to judge the morphology. Individual bursts can contain multiple structures; a burst-level group can hide a short transient superposed on oscillations.

The impulsiveness ranking exposes isolated large transients without a physical template. Start with events **806, 765, 805, and 175**, and inspect saturation-flagged events **497, 784, 1, and 0** separately. Their origin remains unassigned. These inspection examples are not a dust-impact catalog.

## Digital SPAN-electron channel

The user identifies the sixth channel as SPAN-electron counts. The CDF samples are nonnegative integers represented as floats. **All __N__ waveform sums exactly reproduce their stored burst totals.** The samples are treated as increments, not as a cumulative counter. Coded SWEAP status and electron/ion selection masks are preserved without assigning undocumented bit meanings.

There are **__ZEROS__ zero-count** bursts and **__SPARSE__ with 1–31 counts**. These remain browsable but are not assigned a count-shape group. The 32-count cutoff is a pragmatic minimum, not an instrument sensitivity limit. The remaining **__COUNTN__** records split into:

| Group | Description | Records | Core | Medoid |
|---|---|---:|---:|---:|
__COUNTTABLE__

The count partition has mean silhouette **__COUNTSIL__** and median resampling ARI **__COUNTARI__**. Adding count rate changes it appreciably (ARI **__RATEARI__**). Some shape descriptors still depend on counting statistics despite excluding absolute rate from this primary partition.

The viewer also shows a **phase-conditioned expected count trace** and residuals. The expected shape is learned from the other days within the same SWEAP status, with one burst-specific gain to match its approximate total. Residual = (observed − expected)/sqrt(expected + 1). This is an inspection normalization, not a calibrated detector-noise model. Residual partitions are less stable (median resampling ARI **__RESARI__**) and should be regarded as tentative rankings:

| Residual group | Records | Core | Median residual RMS | Medoid |
|---|---:|---:|---:|---:|
__RESTABLE__

The **joint __K__-way partition** uses **__JOINTN__ unsaturated common-configuration records with at least 32 counts**. Its ARI against the analog partition on the same subset is **__JOINTARI__**; analog and raw-count partitions themselves have ARI **__CROSSARI__**. They organize events differently. The catalog therefore retains analog G, raw-count C, residual R, and joint J labels separately; their integer numbering is independent.

## Method and verification

1. **Read-only ingestion.** Each record supplies 32,768 samples at 1.92 MSa/s, or 17.0667 ms of exposure. Arrays are padded to 262,144; only the actual samples enter analysis. All analog-channel existence flags were checked; only channels marked present were read; all loaded waveform values are finite and non-fill. Counts are integer and nonnegative with exact total agreement. TT2000 integers are preserved as decimal strings to avoid browser number-precision loss. Records have distinct file/record provenance and timestamps. Sample timing checks at file boundaries agree with 520/521 ns quantized spacing. Event and CDF record indices are zero-based.
2. **Engineering units.** The frequency-corrected physical-units flag is zero throughout. Electric channels stay in engineering mV and SCM5 in engineering nT, without transfer-function correction. Antenna voltages are not converted to electric fields. All records share 864 kHz low-pass settings and are read without downsampling. Complete alternate records have V5/SCM4 instead of V4/SCM5. Incomplete configurations and missing electron channels are enumerated in [the coverage inventory](coverage.json). These records have no primary G or joint assignment. A separate analysis retains the complete V5/SCM4 configuration; see alternate_configuration_report.md. Their available channels remain browsable. Missing electron channels are labeled absent and excluded from count-reference training and count-based fits; zero-shaped internal preview placeholders are not observed zero counts. Primary canonical feature fields are unavailable, not imputed; actual alternate-channel features and voltage/count associations have separate tables. Thrusting and SCM-calibration flags are zero. The TWTA-on flag is set in __TWTA__ records; it is preserved and included in the adjusted voltage/count comparison. Configuration variables are checked separately from clustering features.
3. **Analog features.** Subtract channel medians. Extract 89 features: normalized crest factor, fourth moment, signed asymmetry, energy duration/asymmetry, RMS-envelope occupancy, energy concentration, zero crossings, autocorrelation, spectral entropy/median/bandwidth/concentration, correlations, envelope correlations, and power-weighted coherence. These features are invariant to positive per-channel gain changes; this was checked numerically on sample records. Five absolute log-RMS amplitudes are saved separately.
4. **Spectral estimation.** Hann-window Welch averages, 4096 samples, 50% overlap, 468.75 Hz bins. Coherence uses the same multi-segment averaging, not a single periodogram. Physical polarization, phase speeds, and wave modes are not inferred from uncorrected channels. Sub-469 Hz spectral detail is poorly resolved by this feature representation, although slow waveform structure remains in temporal descriptors.
5. **Clustering.** Median/IQR scaling with standard-deviation fallback; clip standardized values to ±5. Temporal, spectral, and cross-channel feature families have equal summed squared-distance weight. Ward hierarchical clustering uses the full feature space, not the two-dimensional PCA display. Examine k = 2–10 and select the best silhouette among partitions with at least 20 records per group. The constrained choice is k = __K__; the unconstrained silhouette optimum is k = __BESTK__. This is a declared browsing rule, not proof of a natural class count.
6. **Assignment policy.** Core requires silhouette ≥0.05 and matched-label agreement ≥0.80 across 50 repetitions drawing 80% of records and 80% of features, with scaling refitted. An outlier lies above the empirical 98th percentile of distance to its 15th neighbor: this deliberately flags __OUTLIERS__ records and is not a calibrated anomaly probability. Other unsaturated records remain ambiguous. All __SATFLAGS__ saturation-flagged records are excluded from fitting but fully plotted.
7. **Digital features.** Count dispersion on several aggregation scales, normalized against random rearrangement of empirical single-sample counts; temporal concentration/asymmetry/entropy; spectral shape; and single-sample Fano factor. Forty record/feature resamples assess count groups. Count-to-analog envelope relationships are separate features in the joint analysis.
8. **Cross-channel inspection.** Counts are summed in 128-sample bins (66.67 µs), paired with analog RMS envelopes. Record the largest absolute circular correlation across five channels and ±16 bins. Circular-shift reference p-values and BH-adjusted values are exploratory: scan structure and nonstationarity can invalidate the null. They are not detections of wave–particle interaction and are not physical labels.
9. **Empirical phase baseline.** Divide the observed sample-phase range into 600 bins, use median binned counts from other days with the same SWEAP status, apply a three-bin median smoother and interpolate. No periodic wrap or assumed sweep frequency is imposed. Fit a per-event overall count gain for residual-shape inspection. The reported cross-day R² uses unscaled other-day predictions with matching SWEAP status. Fewer than 20 matching-status training records on other days triggers a flagged cross-status preview, excluded from this R² and residual-count classification. Timing shifts, instrumental scan details, and changing plasma conditions can remain in the residuals.
10. **Line sensitivity.** Detect local spectral excess >8 dB above a 31-bin median in at least 60% of unsaturated records, restrict candidate peaks above 7.5 kHz, and mask ±1 Welch bin around each peak. Recompute four spectral descriptors per channel, retaining temporal/coupling features. This can remove real persistent signals; no time-domain filtering was applied.

## Stability and nuisance checks

ARI = 1 means identical partitions; values near zero indicate little agreement beyond chance. Comparisons keep the same group count to isolate representation changes.

| Perturbation | ARI against analog baseline |
|---|---:|
__ALTTABLE__
| Mask persistent spectral-line bins | __LINEARI__ |

Average linkage produces group sizes __AVGCOUNTS__, collapsing most records into one group. That reinforces the limited separation in this feature representation. Day, onboard selection type, and coded SWEAP status were excluded from the grouping features and cross-tabulated afterward. Group proportions differ by date. This could reflect physical evolution, selection, or instrument-related behavior; no causal interpretation is assigned. Full cross-tabs and association coefficients are in the diagnostic JSON.

Onboard selection limits population interpretation: these results describe the recorded bursts, not unbiased occurrence rates in the solar wind. “Without physics priors” here means no physical labels/templates or event ontology; feature choices, distance metrics, and thresholds still impose statistical assumptions.

## Electron–voltage correlations

A separate [correlation analysis](electron_voltage_report.md) compares voltage RMS with total counts across bursts, and signed voltage means and RMS envelopes with raw and phase-conditioned counts within bursts. Its catalog includes lag searches and exploratory scan-preserving reference comparisons.

## SCM–voltage associations

The [SCM association analysis](scm_voltage_report.md) adds explicit per-event flags and pairwise coefficients for every available voltage with SCM5 or SCM4. Waveform and envelope correlation, coherence outside persistent lines, and line-dominated coherence remain distinct screening categories.

## Coverage and physical identifications

See [the methods and coverage page](methods.html) for E23 date tags, daily file coverage, missing channels, and the PI’s dust-impact reference and high-PC1 regional interpretation. Physical annotations remain separate from unsupervised groups.

## Outputs

Open **atlas.html** in a current browser. It works offline and contains all __N__ events, six-channel previews, filtering, nearest analog neighbors, spectra, count baselines/residuals, and local annotation/export. Full-record displays are 256-bin min/max envelopes, not downsampled waveforms; do not infer frequency from their apparent ripples. Viewer zooms contain 512 actual samples at the strongest normalized electric excursion. Example PNGs use 1024 actual samples and the original Welch frequency bins. Viewer spectra are log-binned previews; classification uses the original bins.

Figure files are in **__FIGURES__**, alongside `data`, as requested. The overview is embedded in the atlas; links to example PNGs point to that local directory. Source CDFs remain unchanged. SHA-256 hashes, features, catalogs, model settings, and scripts accompany the atlas. The HTML contains data-derived previews, not copies of the original CDFs.

The next scientific step is expert inspection of representatives and outliers, plus instrument-informed interpretation of persistent lines and SWEAP scan timing/status. That can guide a second pass while preserving this blind baseline.
