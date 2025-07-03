function patientWalkInfo = ExtractPatientWalkFromFileNames(fileNames, pattern)
    tokens = regexp(fileNames, pattern, 'tokens');
    tokens = vertcat(tokens{:});
    if isempty(tokens)
        patientWalkInfo = [];
        return
    else
        patientWalkInfo = cat(2,cellfun(@(x) str2double(x(1)),tokens),cellfun(@(x) str2double(x(2)),tokens));
    end
end