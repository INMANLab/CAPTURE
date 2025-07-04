# RWAnalysis2 Class Documentation

This document provides an overview of the main methods in the `RWAnalysis2.m` MATLAB class, used for analyzing real-world behavioral and neural data. Each section below links to a specific method with a description and optional usage notes.

---

## getMultData(obj, varargin)

The `getMultData` method processes multi-patient neural recordings and produces two core data structures:

1. **`SegTable`** — a table of 2-second segment-level features used for GLME modeling.
2. **`MultTrans`** — a set of transition-aligned matrices used for event-related spectrogram analysis (e.g. around walk starts/stops).

### 1. SegTable for GLME

Each **row** represents a 2-second segment from a given patient and channel. Each **column** is a feature or label that can be used in a GLME model.

#### Metadata Labels

| Label         | Description |
|--------------|-------------|
| `Pt`          | Patient ID (numeric) |
| `Chan`        | Channel number (1–4) |
| `Walk`        | Walk number or trial index |
| `StopWalk`    | Boolean flag for "stop walk" condition |
| `GoWalk`      | Boolean flag for "go walk" condition |
| `InFlag`      | Inclusion flag (1 = include) |
| `RegionNum`   | Region ID (numeric) |
| `RegionLabel` | Name of brain region (string) |
| `ChanLabel`   | Channel label (string) |

#### Behavioral & Physiological Predictors

| Label     | Description |
|-----------|-------------|
| `Vel`     | Velocity (walking speed), normalized |
| `Fix`     | Mean fixation percentage, normalized |
| `Amb`     | Median ambient light measure, log-normalized |
| `KDE`     | Kernel density estimate, median |
| `EDA`     | Electrodermal activity, normalized |
| `HeadTurn`| Mean of gyroY (head turn), not normalized |

#### Spectral Power Features

Each frequency band appears in two forms:
- `wv*`: Computed using wavelet transforms
- `mt*`: Computed using multitaper spectral estimates

| Label       | Description |
|-------------|-------------|
| `wvTheta`, `mtTheta` | Theta band power (4–8 Hz) |
| `wvAlpha`, `mtAlpha` | Alpha band power (8–12 Hz) |
| `wvBeta`,  `mtBeta`  | Beta band power (12–30 Hz) |
| `wvGamma`, `mtGamma` | Gamma band power (30–70 Hz) |
| `wvHG`,    `mtHG`    | High gamma power (70–120 Hz) |

#### Pepisode (Burst) Features

| Label         | Description |
|---------------|-------------|
| `peTheta`     | Percentage of time with theta bursts |
| `peAlpha`     | Percentage of time with alpha bursts |
| `peBeta`      | Percentage of time with beta bursts |

#### FOOOF-Derived Spectral Features

These features are extracted using the [FOOOF](https://fooof-tools.github.io/fooof/) algorithm, which separates periodic and aperiodic components of the power spectrum.

| Label               | Description |
|---------------------|-------------|
| `peak_freq`         | Frequency of highest spectral peak |
| `peak_height`       | Peak power above aperiodic background |
| `peak_width`        | Width of spectral peak (Hz) |
| `aperiodic_offset`  | Power offset of the aperiodic background |
| `aperiodic_exponent`| Exponent (slope) of the aperiodic component |
| `fit_rsquared`      | R² of FOOOF model fit |
| `fit_error`         | Error of the FOOOF fit |

### 2. MultTrans struct for Spectrogram Analysis

In addition to the segment table, `getMultData` produces event-aligned matrices used for spectrogram visualizations (e.g. via `MultTransSpecGram`).

| Field           | Description |
|------------------|-------------|
| `DT_np`          | Neural power matrix (time × trial × channel) |
| `DT_np_idx`      | Time index of event-aligned windows |
| `DT_xs_idx`      | Velocity-aligned sample index |
| `DT_gz_idx`, `DT_am_idx`, ... | Other behavioral feature-aligned indices |
| `DT_Evnt`        | Event labels (e.g. stop/go) |
| `DT_Desc`        | Event descriptions |
| `DT_StopWalk`, `DT_GoWalk` | Boolean flags for walking transitions |
| `DT_Region`, `DT_Chan`, `DT_Patient` | Metadata arrays for selecting patient/channel subsets for plotting |
| `OL`             | Outlier flags (NaNs, artifactual trials) |
| `OL2`            | Simple IED detection + NaNs combined across all channels/pts |
| `TimeSamp_*`, `TimeSec_*` | Sample and second indices for plotting aligned time axes |
| `FS_*`           | Sampling rates for each data stream |
| `F_wv`, `F_pe`   | Frequency vectors used for wavelet/pepisode transforms |
| `WinSegSec`, `WinTransSec` | Analysis window lengths (in seconds) |

### Summary

| Output         | Description |
|----------------|-------------|
| `obj.MultSeg.SegTable` | Segment-wise feature table used for GLME modeling |
| `obj.MultTrans`        | Event-aligned neural + behavioral data for spectrogram analysis |
| `obj.MultTable`        | Per-patient raw and derived data (used internally) |

This table forms the basis for selecting variables in the GLME GUI. Each variable can be used as a fixed or random effect in the model, and the structure supports multi-level modeling across patients, brain regions, and channels.
