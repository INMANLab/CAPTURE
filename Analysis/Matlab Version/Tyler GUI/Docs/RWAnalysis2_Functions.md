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

| **Field**         | **Description** |
|------------------|-----------------|
| `DT_np_idx`       | Index of each segment's center sample (for neural power) |
| `DT_xs_idx`       | Sample index for velocity-aligned events |
| `DT_xs_vchg`      | Percent change in velocity for each trial |
| `DT_gz_idx`       | Sample index aligned to fixation events |
| `DT_am_idx`       | Sample index aligned to ambient condition changes |
| `DT_kd_idx`       | Sample index for kernel density estimate (position) |
| `DT_wv_idx`       | Sample index aligned to wavelet power data |
| `DT_bi_idx`       | Sample index for biometric (e.g., EDA) data |
| `DT_pe_idx`       | Sample index for pepisode (burst detection) data |
| `DT_im_idx`       | Sample index for inertial movement (IMU/gyro) data |
| `Walk`            | Walk number associated with each transition trial |
| `NSamp`           | Number of samples in each transition window |
| `Evnt`            | Event label (e.g., stop or go) |
| `Desc`            | Event description string |
| `StopWalk`        | Boolean flag indicating "stop walk" transitions |
| `GoWalk`          | Boolean flag indicating "go walk" transitions |
| `OL`              | Outlier flags: NaN or -500 (used to exclude noisy trials) |
| `OL2`             | Combined artifact detection: NaNs and z-scored outlier power (IED detection) |
| `Region`          | Brain region numeric ID |
| `RegionLabel`     | Brain region name string |
| `Chan`            | Channel ID for each trial |
| `ChanLabel`       | Channel label string (e.g., `AntHipp-a`) |
| `Patient`         | Patient ID associated with each transition |
| `TimeSamp_np`     | Sample indices for neural power (wavelet) |
| `TimeSec_np`      | Same as above, but in seconds |
| `TimeSamp_xs`     | Sample indices for velocity time series |
| `TimeSec_xs`      | Velocity time axis in seconds |
| `TimeSamp_gz`     | Sample indices for fixation time series |
| `TimeSec_gz`      | Fixation time axis in seconds |
| `TimeSamp_am`     | Sample indices for ambient signal |
| `TimeSec_am`      | Ambient signal time axis in seconds |
| `TimeSamp_kd`     | Sample indices for KDE signal |
| `TimeSec_kd`      | KDE time axis in seconds |
| `TimeSamp_wv`     | Sample indices for wavelet power |
| `TimeSec_wv`      | Wavelet time axis in seconds |
| `TimeSamp_bi`     | Sample indices for biometric signal |
| `TimeSec_bi`      | Biometric time axis in seconds |
| `TimeSamp_pe`     | Sample indices for pepisode signal |
| `TimeSec_pe`      | Pepisode time axis in seconds |
| `TimeSamp_im`     | Sample indices for inertial signal |
| `TimeSec_im`      | Inertial/IMU time axis in seconds |
| `FS_np`           | Sampling rate for neural power (wavelet) |
| `FS_xs`           | Sampling rate for velocity data |
| `FS_gz`           | Sampling rate for fixation data |
| `FS_am`           | Sampling rate for ambient signal |
| `FS_kd`           | Sampling rate for KDE signal |
| `FS_wv`           | Sampling rate for wavelet power |
| `FS_bi`           | Sampling rate for biometric data (e.g., EDA) |
| `FS_pe`           | Sampling rate for pepisode burst detection |
| `FS_im`           | Sampling rate for IMU/inertial data |
| `F_bi`            | Field labels for biometric channels |
| `F_im`            | Field labels for inertial movement channels |
| `F_pe`            | Frequency bands used in pepisode analysis |
| `Freq`            | Frequency vector for wavelet transform |
| `WinSegSec`       | Window size in seconds for segment analysis and smoothing |
| `WinTransSec`     | Window size in seconds for transition-aligned analysis |

### Summary

| Output         | Description |
|----------------|-------------|
| `obj.MultSeg.SegTable` | Segment-wise feature table used for GLME modeling |
| `obj.MultTrans`        | Event-aligned neural + behavioral data for spectrogram analysis |
| `obj.MultTable`        | Per-patient raw and derived data (used internally) |

This table forms the basis for selecting variables in the GLME GUI. Each variable can be used as a fixed or random effect in the model, and the structure supports multi-level modeling across patients, brain regions, and channels.


---
## saveMultData(obj,varargin)
Saves the processed data structures to a `.mat` file for persistent storage.

#### Save format:

| Variable       | Description |
|----------------|-------------|
| `multseg`      | The full segment-level table and related metadata (`obj.MultSeg`) used for GLME modeling |
| `multtrans`    | Transition-aligned neural and behavioral data (`obj.MultTrans`) used for spectrograms |
| `multtable`    | Per-patient processed data (`obj.MultTable`) used internally or for debug |

```matlab
save(obj.AnalysisFile,'multseg','multtrans','multtable','-v7.3');
```

---
## loadMultData(obj,varargin)

Loads The `loadMultData` function restores previously saved analysis data from disk if the .mat file exists, allowing you to resume work without re-running the full `getMultData` pipeline.

| Saved Field   | Loaded Into        | Description |
|---------------|--------------------|-------------|
| `multseg`     | `obj.MultSeg`      | Segment-level table and metadata for GLME modeling |
| `multtrans`   | `obj.MultTrans`    | Transition-aligned data for spectrogram analysis |
| `multtable`   | `obj.MultTable`    | Per-patient data tables used internally |

