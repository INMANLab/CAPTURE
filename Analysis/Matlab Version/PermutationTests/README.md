# MATLAB Permutation Tests (`+permutation`)

This folder contains a modular MATLAB package for participant-level permutation tests using trial-level data.

## What This Implements

- **One-sample test vs zero**  
  Test participant-mean values against zero, with trial-level sign-flip permutations within each participant.
- **Two-sample test**  
  Test paired trial-level conditions (`group1` vs `group2`) by shuffling condition labels within each participant.
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

All statistical comparisons are computed at the **participant level** (average within participant, then aggregate across participants).

---

## Public API

### 1) One-sample vs zero

```matlab
[nullDist, p, stat, info] = permutation.OneSampleVsZero( ...
    group1, patientList, nPerm=5000, Seed=42);
```

Returns:
- `nullDist`: null distribution (`nPerm x 1`)
- `p`: two-sided Monte Carlo p-value
- `stat`: observed statistic
- `info`: struct with participant means and metadata (`nPerm`, `participantIds`, feasibility info)

### 2) Two-sample

```matlab
[nullDist, p, stat, info] = permutation.TwoSample( ...
    group1, group2, patientList, nPerm=5000, Seed=42);
```

Returns same output structure as one-sample.

### 3) Plot null distribution

```matlab
permutation.PlotNullDistribution(nullDist, stat, PValue=p);
```

### 4) Plot participant averages (boxplot + overlay summary)

```matlab
permutation.PlotParticipantAveragesBoxplot(p, stat, info);
```

This function reads participant averages from `info` and overlays:
- p-value
- number of permutations (`info.nPerm`)
- statistic value

### 5) Feasibility utilities

Detailed log-scale feasibility:

```matlab
out = permutation.Feasibility('twoSample', patientList, nPerm, group1, group2);
disp(out.message);
```

Quick participant-wise risk check:

```matlab
[riskFlag, validationRisk] = permutation.FeasibilityCheck('twoSample', patientList, nPerm);
```

---

## Example Script

See `PermutationTestScript.m` for an end-to-end runnable example using `RWA_PowerBox_Example.mat`.

Run:

```matlab
PermutationTestScript
```

---

## Notes

- Ensure `Functions` is on your MATLAB path so package calls like `permutation.OneSampleVsZero(...)` resolve.
- Two-sample test assumes aligned trial rows between `group1` and `group2` for each participant.
