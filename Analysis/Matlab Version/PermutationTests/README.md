# MATLAB Permutation Tests (`+permutation`)

This folder contains a modular MATLAB package for participant-level permutation tests using trial-level data, with an optional **(participant × channel)** nesting layer.

## What This Implements

- **One-sample test vs zero**  
  Test aggregated values against zero, with trial-level sign-flip permutations within each participant (or within each participant–channel stratum when `patientChannel` is supplied).
- **Two-sample test**  
  Test paired trial-level conditions (`group1` vs `group2`) by shuffling condition labels within each participant (or within each stratum when channels are supplied).
- **Null distribution plotting**
- **Participant-average boxplot plotting with test summary overlay**
- **Permutation-space feasibility checks**

Core package lives in `+permutation`.

---

## Data Convention

Inputs are trial-level vectors:

- `group1` : `N x 1` numeric vector
- `group2` : `N x 1` numeric vector (two-sample only)
- `patientList` : `N x 1` participant IDs matching each trial row
- `nPerm` : number of permutations (default `1000`)
- **`patientChannel`** (optional) : `N x 1` channel labels aligned with each trial. When provided, each **stratum** is a unique `(participant, channel)` pair; sign-flips / label shuffles occur **only within that stratum**.

### Aggregation when using channels

Stratum means \(\mu_{p,c}\) are computed per `(participant, channel)`. The scalar test statistic is then:

- **`Aggregation="participant"`** (default): `mean_p( mean_c( μ_{p,c} | same p ) )` — average across channels within each participant, then across participants.
- **`Aggregation="channel"`**: `mean_c( mean_p( μ_{p,c} | same c ) )` — average across participants within each channel, then across channels.

For **two-sample**, the same aggregation is applied to **stratum-wise** differences \(\mu^{(1)}_{p,c} - \mu^{(2)}_{p,c}\).

If `patientChannel` is omitted, behavior matches the original participant-only pipeline (aggregation name-value is ignored).

---

## Public API

### 1) One-sample vs zero

```matlab
[nullDist, p, stat, info] = permutation.OneSampleVsZero( ...
    group1, patientList, nPerm=5000, Seed=42);
```

With nested strata:

```matlab
[nullDist, p, stat, info] = permutation.OneSampleVsZero( ...
    group1, patientList, nPerm=5000, Seed=42, ...
    patientChannel=patientChannel, Aggregation="participant");
```

Returns:

- `nullDist`: null distribution (`nPerm x 1`)
- `p`: two-sided Monte Carlo p-value
- `stat`: observed statistic
- `info`: struct with means, `nPerm`, `participantIds`, feasibility; if nested, also `stratumMeansObserved`, `stratumParticipantIds`, `stratumChannelIds`, `Aggregation`

### 2) Two-sample

```matlab
[nullDist, p, stat, info] = permutation.TwoSample( ...
    group1, group2, patientList, nPerm=5000, Seed=42);
```

Nested:

```matlab
[nullDist, p, stat, info] = permutation.TwoSample( ...
    group1, group2, patientList, nPerm=5000, Seed=42, ...
    patientChannel=patientChannel, Aggregation="channel");
```

### 3) Plot null distribution

```matlab
permutation.PlotNullDistribution(nullDist, stat, PValue=p);
```

### 4) Plot participant averages (boxplot + overlay summary)

```matlab
permutation.PlotParticipantAveragesBoxplot(p, stat, info);
```

Optional nested overlay: one point per `(participant, channel)` stratum mean, with **one marker shape per participant** (e.g. two `^` for two channels, four `s` for four channels). X-position is spread by channel within each box group.

```matlab
permutation.PlotParticipantAveragesBoxplot(p, stat, info, ...
    PlotPerChannelStratumPoints=true, ShowParticipantLegend=true);
% Optional: only triangles and squares in cycle — ParticipantMarkerCycle=["^","s"]
```

This function reads participant-level values from `info` (for nested runs: mean of stratum means across channels per participant, for visualization) and overlays p-value, `nPerm`, statistic, and nested/aggregation notes when present.

### 5) Feasibility utilities

Detailed log-scale feasibility:

```matlab
out = permutation.Feasibility('twoSample', patientList, nPerm, group1, group2);
disp(out.message);

% Optional 6th argument: trial-level channel labels (same length as patientList)
out = permutation.Feasibility('twoSample', patientList, nPerm, group1, group2, patientChannel);
```

One-sample (fourth argument is `group1`; unused `group2` defaults empty):

```matlab
out = permutation.Feasibility('oneSample', patientList, nPerm, group1);
out = permutation.Feasibility('oneSample', patientList, nPerm, group1, [], patientChannel);
```

Quick per-stratum risk check:

```matlab
[riskFlag, validationRisk] = permutation.FeasibilityCheck('twoSample', patientList, nPerm);
[riskFlag, validationRisk] = permutation.FeasibilityCheck('twoSample', patientList, nPerm, patientChannel);
```

---

## Example Script

See `PermutationTestScript.m` for an end-to-end runnable example using `RWA_PowerBox_Example.mat` (and optional nested mode using `channels` from the file or a synthetic `patientChannel`).

Run:

```matlab
PermutationTestScript
```

---

## Notes

- Add the folder that contains `+permutation` to your MATLAB path so calls like `permutation.OneSampleVsZero(...)` resolve.
- Two-sample test assumes aligned trial rows between `group1` and `group2` for each participant (and within each channel when `patientChannel` is used).
