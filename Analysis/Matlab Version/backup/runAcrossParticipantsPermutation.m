function result = runAcrossParticipantsPermutation(group1, group2, patientList, testType, nPerm, statType, tail, doFDR, q)
%RUNACROSSPARTICIPANTSPERMUTATION Cohort test with trial-within-participant permutations.
%
%   This entry point now delegates to the traceable cohort pipeline:
%   oneSampleVsZero_trialPermutedCohortMean or twoGroupPaired_trialPermutedCohortMean
%   on the full trial-level vectors (not pre-aggregated participant means).
%
%   statType is accepted for API compatibility but ignored (cohort statistic is
%   always mean across participants of participant summaries).

    if nargin < 4
        error('permutation:runAcrossParticipantsPermutation:MissingInputs', ...
            'Required inputs: group1, group2, patientList, testType.');
    end
    if nargin < 5 || isempty(nPerm), nPerm = 1000; end
    if nargin < 6, statType = []; end %#ok<NASGU>
    if nargin < 7 || isempty(tail), tail = 'two-sided'; end
    if nargin < 8 || isempty(doFDR), doFDR = false; end
    if nargin < 9 || isempty(q), q = 0.05; end

    testType = validatestring(lower(char(testType)), {'one-sample-vs-zero','two-sample-paired'});
    tail = validatestring(lower(char(tail)), {'two-sided','right','left'});
    doFDR = logical(doFDR);

    switch testType
        case 'one-sample-vs-zero'
            overall = permutation.oneSampleVsZero_trialPermutedCohortMean( ...
                group1, patientList, nPerm, tail);
        case 'two-sample-paired'
            overall = permutation.twoGroupPaired_trialPermutedCohortMean( ...
                group1, group2, patientList, nPerm, tail);
    end

    if doFDR
        [pAdj, sigMask, criticalP] = permutation.applyFDR_BH(overall.pValue, q);
        overall.pValueAdj = pAdj;
        overall.isSignificant = sigMask;
        overall.criticalP = criticalP;
    else
        overall.pValueAdj = NaN;
        overall.isSignificant = false;
        overall.criticalP = NaN;
    end

    aggG1 = permutation.aggregateParticipantMeans(group1, patientList);
    summaryG2 = [];
    if strcmp(testType, 'two-sample-paired')
        aggG2 = permutation.aggregateParticipantMeans(group2, patientList);
        summaryG2 = aggG2.participantMeans;
    end

    overall.metadata.summaryRule = 'trialPerm_withinParticipant_cohortMean';

    result = struct( ...
        'testType', testType, ...
        'nPerm', nPerm, ...
        'statType', 'cohortMeanOfParticipantSummaries', ...
        'tail', tail, ...
        'patientSummary', struct('group1', aggG1.participantMeans, 'group2', summaryG2), ...
        'overall', overall);
end
