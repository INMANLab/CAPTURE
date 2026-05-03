# Trial-level permutation, cohort statistic (MATLAB)

This package implements **cohort-level permutation tests** where:

1. You have **trial-level** observations (length `N`) with a **participant label** per row (`patientList`, `N×1`).
2. The **test statistic** is computed from **one summary scalar per participant** (for example mean of trials within that participant, or mean of paired differences within that participant).
3. The **cohort statistic** is the **mean across participants** of those per-participant summaries.
4. The **null distribution** is built by **permuting at trial level within each participant only** (independently across participants), then recomputing participant summaries and the cohort mean, for each of `nPerm` replicates.

This matches the design: *trial-level randomization within participant → participant-level averages → one cohort value per permutation*.

---

## Primary entry points (recommended)

| Test | Function |
|------|----------|
| One group vs zero | `permutation.oneSampleVsZero_trialPermutedCohortMean(group1, patientList, nPerm, tail)` |
| Two groups (paired trials within participant) | `permutation.twoGroupPaired_trialPermutedCohortMean(group1, group2, patientList, nPerm, tail)` |

Stable aliases:

- `permutation.oneSampleVsZero_test(...)` → same as `oneSampleVsZero_trialPermutedCohortMean`
- `permutation.twoSamplePaired_test(...)` → same as `twoGroupPaired_trialPermutedCohortMean`

---

## Exact definitions

### One-sample vs zero

For each participant `p` (stable order `unique(patientList,'stable')`):

\[
m_p = \mathrm{mean}(\text{group1 trials of } p,\ \text{omitnan})
\]

Cohort statistic:

\[
T = \mathrm{mean}_p(m_p,\ \text{omitnan})
\]

**Null replicate `b`:** draw independent trial-wise sign flips `±1` within each participant, multiply trial values, recompute all `m_p^{(b)}`, then \(T_b = \mathrm{mean}_p(m_p^{(b)})\).

### Two-group paired (row-matched trials)

For each participant `p`:

\[
d_p = \mathrm{mean}(\text{group1}_p - \text{group2}_p,\ \text{omitnan})
\]

Cohort statistic:

\[
T = \mathrm{mean}_p(d_p,\ \text{omitnan})
\]

**Null replicate `b`:** within each participant, independently at each trial either **keep** or **swap** the `(group1, group2)` labels; recompute \(d_p^{(b)}\), then \(T_b = \mathrm{mean}_p(d_p^{(b)})\).

### P-value

See `permutationPValue.m`. Default two-sided:

\[
p = \frac{1 + \sum_{b=1}^{B} \mathbf{1}(|T_b| \ge |T_{\text{obs}}|)}{B + 1}
\]

---

## Traceable building blocks

| File | Role |
|------|------|
| [`+permutation/generateTrialLevelPermutationLabels.m`](+permutation/generateTrialLevelPermutationLabels.m) | Builds per-participant `nTrials × nPerm` label matrices (column 1 = identity). |
| [`+permutation/applyTrialLevelPermutation.m`](+permutation/applyTrialLevelPermutation.m) | Applies one permutation column to full `N×1` vectors. |
| [`+permutation/aggregateParticipantMeans.m`](+permutation/aggregateParticipantMeans.m) | Per-participant mean of `group1`. |
| [`+permutation/aggregateParticipantPairedMeanDiff.m`](+permutation/aggregateParticipantPairedMeanDiff.m) | Per-participant mean of `group1 - group2`. |
| [`+permutation/cohortMeanStatistic.m`](+permutation/cohortMeanStatistic.m) | Mean across the `P×1` participant vector. |
| [`+permutation/permutationPValue.m`](+permutation/permutationPValue.m) | Shared p-value counting rule. |

---

## Other utilities

- **FDR (Benjamini–Hochberg):** [`+permutation/applyFDR_BH.m`](+permutation/applyFDR_BH.m) — optional; not part of the single cohort statistic definition.
- **Plot cohort null + summaries:** [`+permutation/plotPermutationResult.m`](+permutation/plotPermutationResult.m) — pass `patientList` as 5th argument to plot **participant-level** summaries on the right panel.

---

## Quick example

```matlab
rng(1);
patientList = repelem((1:5)', 20);   % 5 participants, 20 trials each
group1 = 0.2 + randn(numel(patientList), 1);
group2 = randn(numel(patientList), 1);

res1 = permutation.oneSampleVsZero_trialPermutedCohortMean(group1, patientList, 2000, 'two-sided');
res2 = permutation.twoGroupPaired_trialPermutedCohortMean(group1, group2, patientList, 2000, 'two-sided');

disp(res1.obsStat), disp(res1.pValue);
disp(res2.obsStat), disp(res2.pValue);
```

Run the bundled demo:

```matlab
PermutationTestScript
```

---

## Inputs

- `group1`, `group2`: `N×1` double. For two-group tests, same length; trials **row-matched** within participant.
- `patientList`: `N×1` participant ID per row (numeric, categorical, string, or char column — column shape required).
- `nPerm`: default `1000` in functions that default it.
- `tail`: `'two-sided'`, `'right'`, or `'left'`.

---

## Outputs (result struct)

Core fields:

- `nullDist` (`nPerm×1`)
- `obsStat` (scalar)
- `pValue`
- `tail`, `nPerm`, `testType`, `statType`
- `labelStruct` (full label object for inspection)
- `metadata` (e.g. `participantMeansObserved` / `participantDiffObserved`, `patientIDs`, `formula`)
