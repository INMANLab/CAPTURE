function result = runAcrossParticipantsPermutation(group1, group2, patientList, testType, nPerm, statType, tail, doFDR, q)
%RUNACROSSPARTICIPANTSPERMUTATION Run permutation tests at participant level.
%   result = permutation.runAcrossParticipantsPermutation(...)
%   first summarizes each participant to one value and then runs one
%   permutation test across participant summaries.

    if nargin < 4
        error('permutation:runAcrossParticipantsPermutation:MissingInputs', ...
            'Required inputs: group1, group2, patientList, testType.');
    end
    if nargin < 5 || isempty(nPerm), nPerm = 1000; end
    if nargin < 6 || isempty(statType), statType = 'mean'; end
    if nargin < 7 || isempty(tail), tail = 'two-sided'; end
    if nargin < 8 || isempty(doFDR), doFDR = true; end
    if nargin < 9 || isempty(q), q = 0.05; end

    validateattributes(group1, {'numeric'}, {'column'}, mfilename, 'group1');
    if ~isempty(group2)
        validateattributes(group2, {'numeric'}, {'column'}, mfilename, 'group2');
    end
    validateattributes(nPerm, {'numeric'}, {'scalar','integer','>=',1}, mfilename, 'nPerm');
    validateattributes(q, {'numeric'}, {'scalar','>',0,'<',1}, mfilename, 'q');
    doFDR = logical(doFDR);

    testType = validatestring(lower(char(testType)), {'one-sample-vs-zero','two-sample-paired'});
    statType = validatestring(lower(char(statType)), {'mean','median','t'});
    tail = validatestring(lower(char(tail)), {'two-sided','right','left'});

    if numel(group1) ~= numel(patientList)
        error('permutation:runAcrossParticipantsPermutation:DimensionMismatch', ...
            'group1 and patientList must have the same length.');
    end

    if strcmpi(testType, 'two-sample-paired')
        if isempty(group2) || numel(group2) ~= numel(group1)
            error('permutation:runAcrossParticipantsPermutation:PairedInputError', ...
                'For paired tests, group2 must be provided and match group1 length.');
        end
    end

    [uniquePatients, ~, idxPatient] = unique(patientList, 'stable');
    nPatients = numel(uniquePatients);

    summaryG1 = nan(nPatients, 1);
    summaryG2 = nan(nPatients, 1);

    for iPatient = 1:nPatients
        mask = idxPatient == iPatient;
        summaryG1(iPatient) = mean(group1(mask), 'omitnan');
        if strcmpi(testType, 'two-sample-paired')
            summaryG2(iPatient) = mean(group2(mask), 'omitnan');
        end
    end

    switch testType
        case 'one-sample-vs-zero'
            overall = permutation.oneSampleVsZero_test(summaryG1, nPerm, statType, tail);
        case 'two-sample-paired'
            overall = permutation.twoSamplePaired_test(summaryG1, summaryG2, nPerm, statType, tail);
    end

    % Optional BH on a single test is identity-like; kept for API symmetry.
    if doFDR
        [pAdj, sigMask, criticalP] = permutation.applyFDR_BH(overall.pValue, q);
        overall.pValueAdj = pAdj;
        overall.isSignificant = sigMask;
        overall.criticalP = criticalP;
    end

    overall.metadata.patientIDs = uniquePatients;
    overall.metadata.nPatients = nPatients;
    overall.metadata.summaryRule = 'mean-within-patient';

    result = struct( ...
        'testType', testType, ...
        'nPerm', nPerm, ...
        'statType', statType, ...
        'tail', tail, ...
        'patientSummary', struct('group1', summaryG1, 'group2', summaryG2), ...
        'overall', overall);
end
