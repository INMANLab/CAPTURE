# Modular MATLAB Permutation Tests

This folder contains a modular MATLAB implementation for permutation testing with:

- one-group vs zero tests
- two-group paired tests
- participant-level and across-participant workflows
- Benjamini-Hochberg false discovery rate (FDR) control for multiple tests

All main functions are in the package folder `+permutation`.

---

## Contents

- [`+permutation/oneSampleVsZero_test.m`](+permutation/oneSampleVsZero_test.m)
- [`+permutation/twoSamplePaired_test.m`](+permutation/twoSamplePaired_test.m)
- [`+permutation/runParticipantPermutation.m`](+permutation/runParticipantPermutation.m)
- [`+permutation/runAcrossParticipantsPermutation.m`](+permutation/runAcrossParticipantsPermutation.m)
- [`+permutation/OneSample_ParticipantLevel.m`](+permutation/OneSample_ParticipantLevel.m)
- [`+permutation/OneSample_AcrossParticipants.m`](+permutation/OneSample_AcrossParticipants.m)
- [`+permutation/generatePermutationLabels.m`](+permutation/generatePermutationLabels.m)
- [`+permutation/computeStatistic.m`](+permutation/computeStatistic.m)
- [`+permutation/applyFDR_BH.m`](+permutation/applyFDR_BH.m)
- [`+permutation/plotPermutationResult.m`](+permutation/plotPermutationResult.m)
- [`+permutation/plotParticipantPermutationSummary.m`](+permutation/plotParticipantPermutationSummary.m)
- Backward-compatible wrappers:
  - [`+permutation/Twosample_test.m`](+permutation/Twosample_test.m)
  - [`+permutation/TwoSample_ParticipantLevel.m`](+permutation/TwoSample_ParticipantLevel.m)
- Example script:
  - [`PermutationTestScript.m`](PermutationTestScript.m)

---

## Quick Start

From MATLAB, set current folder to this directory and run:

```matlab
PermutationTestScript
```

This script runs deterministic smoke tests on synthetic data and an optional legacy dataset demo.

---

## Data Format

Expected column vectors:

- `group1`: `N x 1` numeric observations
- `group2`: `N x 1` numeric observations (required for paired two-sample tests)
- `patientList`: `N x 1` participant identifiers aligned with observations
- `nPerm`: number of permutations (default: `1000`)

NaN handling:

- Single-test functions remove invalid samples (`omitnan` behavior)
- Across-participant summaries use participant means with `omitnan`

---

## Core API

### 1) One-sample vs zero

```matlab
result = permutation.oneSampleVsZero_test(group1, nPerm, statType, tail)
```

Use sign-flip permutations against zero.

### 2) Paired two-sample test

```matlab
result = permutation.twoSamplePaired_test(group1, group2, nPerm, statType, tail)
```

Use within-pair random swaps (equivalent to sign flips on paired differences).

### 3) Per-participant tests

```matlab
result = permutation.runParticipantPermutation( ...
    group1, group2, patientList, testType, nPerm, statType, tail, doFDR, q)
```

Splits data by participant and runs one test per participant.

### 4) Across-participant test

```matlab
result = permutation.runAcrossParticipantsPermutation( ...
    group1, group2, patientList, testType, nPerm, statType, tail, doFDR, q)
```

Builds one summary value per participant, then runs a cohort-level permutation test.

### 5) FDR correction utility

```matlab
[pAdj, isSignificant, criticalP] = permutation.applyFDR_BH(pValues, q)
```

Applies Benjamini-Hochberg correction.

### 6) One-sample convenience entry points

```matlab
% one-sample-vs-zero within each participant
resWithin = permutation.OneSample_ParticipantLevel( ...
    patientList, group1, nPerm, statType, tail, true, 0.05);

% one-sample-vs-zero across participants
resAcross = permutation.OneSample_AcrossParticipants( ...
    patientList, group1, nPerm, statType, tail, true, 0.05);
```

These wrappers avoid passing an unused `group2`.

### 7) Plotting utility

```matlab
fig = permutation.plotPermutationResult(result, group1, group2, plotTitle)
```

Features:

- plots `nullDist` histogram and superimposes observed statistic (`obsStat`)
- plots boxplot of observed values (`group1`, or `group1/group2` for paired)
- annotates test metadata and p-values directly on the figure

### 8) Participant-level plotting utility

```matlab
fig = permutation.plotParticipantPermutationSummary( ...
    participantResult, group1, group2, patientList, plotTitle)
```

Features:

- left panel: participant-wise boxplot-style distributions
- right panel: per-participant raw p-values, adjusted p-values, alpha line
- highlights FDR-significant participants

---

## Parameters

- `testType`:
  - `'one-sample-vs-zero'`
  - `'two-sample-paired'`
- `statType`:
  - `'mean'`
  - `'median'`
  - `'t'`
- `tail`:
  - `'two-sided'`
  - `'right'`
  - `'left'`
- `doFDR`:
  - logical (`true`/`false`)
- `q`:
  - FDR target level, default `0.05`

---

## Standard Output Contract

Single-test functions return a struct with:

- `nullDist`: permutation null distribution
- `obsStat`: observed statistic
- `pValue`: raw permutation p-value
- `pValueAdj`: adjusted p-value (filled when FDR applied)
- `isSignificant`: significance after FDR (when applied)
- `criticalP`: BH critical p threshold (when applicable)
- `statType`, `tail`, `nPerm`, `testType`
- `metadata`: additional context (`nObs`, participant info, etc.)

Orchestrators include additional aggregation fields:

- `runParticipantPermutation`: `patientIDs`, `perPatient`
- `runAcrossParticipantsPermutation`: `patientSummary`, `overall`

---

## Examples

### Example A: One-sample vs zero

```matlab
rng(7);
group1 = 0.3 + randn(120,1);
res = permutation.oneSampleVsZero_test(group1, 2000, 'mean', 'two-sided');
disp(res.pValue)
```

### Example B: Paired two-sample

```matlab
rng(7);
group1 = 0.5 + randn(120,1);
group2 = randn(120,1);
res = permutation.twoSamplePaired_test(group1, group2, 2000, 't', 'two-sided');
disp([res.obsStat, res.pValue])
```

### Example C: Per-participant with FDR

```matlab
rng(7);
nPatients = 8;
nObsEach = 15;
N = nPatients * nObsEach;
patientList = repelem((1:nPatients)', nObsEach);
group1 = 0.4 + randn(N,1);
group2 = randn(N,1);

res = permutation.runParticipantPermutation( ...
    group1, group2, patientList, ...
    'two-sample-paired', 1000, 'mean', 'two-sided', true, 0.05);

pRaw = arrayfun(@(s) s.pValue, res.perPatient);
pAdj = arrayfun(@(s) s.pValueAdj, res.perPatient);
disp(table((1:numel(pRaw))', pRaw(:), pAdj(:), ...
    'VariableNames', {'PatientIdx','PValueRaw','PValueFDR'}))
```

### Example D: Across-participant aggregate test

```matlab
res = permutation.runAcrossParticipantsPermutation( ...
    group1, group2, patientList, ...
    'two-sample-paired', 2000, 'mean', 'two-sided', true, 0.05);
disp(res.overall.pValue)
```

---

## Backward Compatibility

Legacy names are preserved as wrappers:

- `permutation.Twosample_test(...)` -> `permutation.twoSamplePaired_test(...)`
- `permutation.TwoSample_ParticipantLevel(...)` -> `permutation.runParticipantPermutation(...)`

---

## Notes and Best Practices

- Use larger `nPerm` (for example `5000` to `10000`) for stable tail p-values.
- For paired designs, ensure `group1`, `group2`, and `patientList` are aligned row-wise.
- If you run many tests (participants/features), keep `doFDR = true`.
- Set random seed (`rng`) to make results reproducible across runs.

---

## Troubleshooting

- **Dimension mismatch errors**: check vector lengths and column shape.
- **No valid data errors**: verify NaN prevalence after filtering.
- **Unexpected p-values**: confirm test direction (`tail`) and statistic (`statType`).
- **FDR fields unchanged**: FDR only changes values meaningfully when multiple tests are corrected together.
