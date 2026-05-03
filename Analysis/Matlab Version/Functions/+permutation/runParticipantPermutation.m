function result = runParticipantPermutation(group1, group2, patientList, testType, nPerm, statType, tail, doFDR, q)
%RUNPARTICIPANTPERMUTATION Run permutation tests separately per participant.
%   result = permutation.runParticipantPermutation(...)
%   splits data by patientList and runs either one-sample-vs-zero or
%   paired two-sample tests for each participant.

    if nargin < 4
        error('permutation:runParticipantPermutation:MissingInputs', ...
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
        error('permutation:runParticipantPermutation:DimensionMismatch', ...
            'group1 and patientList must have the same length.');
    end

    if strcmpi(testType, 'two-sample-paired')
        if isempty(group2) || numel(group2) ~= numel(group1)
            error('permutation:runParticipantPermutation:PairedInputError', ...
                'For paired tests, group2 must be provided and match group1 length.');
        end
    end

    [uniquePatients, ~, idxPatient] = unique(patientList, 'stable');
    nPatients = numel(uniquePatients);
    perPatient = cell(nPatients, 1);
    pVals = nan(nPatients, 1);

    for iPatient = 1:nPatients
        mask = idxPatient == iPatient;
        currentG1 = group1(mask);

        switch testType
            case 'one-sample-vs-zero'
                perPatient{iPatient} = permutation.oneSampleVsZero_test(currentG1, nPerm, statType, tail);
            case 'two-sample-paired'
                currentG2 = group2(mask);
                perPatient{iPatient} = permutation.twoSamplePaired_test(currentG1, currentG2, nPerm, statType, tail);
        end

        perPatient{iPatient}.metadata.patientID = uniquePatients(iPatient);
        pVals(iPatient) = perPatient{iPatient}.pValue;
    end

    perPatient = vertcat(perPatient{:});

    if doFDR && nPatients > 1
        [pAdj, sigMask, criticalP] = permutation.applyFDR_BH(pVals, q);
        for iPatient = 1:nPatients
            perPatient(iPatient).pValueAdj = pAdj(iPatient);
            perPatient(iPatient).isSignificant = sigMask(iPatient);
            perPatient(iPatient).criticalP = criticalP;
        end
    end

    result = struct( ...
        'testType', testType, ...
        'nPerm', nPerm, ...
        'statType', statType, ...
        'tail', tail, ...
        'patientIDs', uniquePatients, ...
        'perPatient', perPatient);
end
