function [group1Perm, group2Perm] = permuteLabelsWithinParticipants( ...
    mode, group1, patientList, group2, rngState, patientChannel)
%PERMUTELABELSWITHINPARTICIPANTS One trial-level null draw (within participant).
%
%   group1Perm = permuteLabelsWithinParticipants("oneSample", group1, patientList)
%   [group1Perm, group2Perm] = permuteLabelsWithinParticipants( ...
%       "twoSample", group1, patientList, group2)
%
%   twoSample : pool [group1; group2] trials per participant (paired rows),
%               random partition of size n into pseudo-group1 and pseudo-group2.
%   oneSample : multiply each trial by i.i.d. Rademacher (+1/-1) within participant.
%
%   rngState : optional RandStream or [] to use default global stream.
%   patientChannel : optional Nx1; if nonempty, permute within (participant, channel).

arguments
    mode
    group1 (:,1) {mustBeNumeric}
    patientList (:,1)
    group2 {mustBeNumeric} = zeros(0, 1)
    rngState = []
    patientChannel (:,1) = []
end

mode = validatestring(mode, {'twoSample', 'oneSample'});

if isempty(rngState)
    rngStream = RandStream.getGlobalStream();
else
    rngStream = rngState;
end

if isempty(patientChannel)
    [~, indexCell] = permutation.participantTrialIndices(patientList);
else
    if numel(patientChannel) ~= numel(patientList)
        error("permuteLabelsWithinParticipants:patientChannel must match patientList length.");
    end
    [~, ~, indexCell] = permutation.participantChannelStrata(patientList, patientChannel);
end
U = numel(indexCell);

if strcmp(mode, 'twoSample')
    if isempty(group2)
        error("permuteLabelsWithinParticipants:twoSample requires group2.");
    end
    group2 = group2(:);
    if ~isequal(size(group1), size(group2))
        error("permuteLabelsWithinParticipants:group1 and group2 must match in size.");
    end
    group1Perm = zeros(size(group1));
    group2Perm = zeros(size(group2));
    for k = 1:U
        idx = indexCell{k};
        nPair = numel(idx);
        if nPair == 0
            continue
        end
        pool = [group1(idx); group2(idx)];
        ord = randperm(rngStream, 2 * nPair);
        pool = pool(ord);
        group1Perm(idx) = pool(1:nPair);
        group2Perm(idx) = pool(nPair + 1:end);
    end
else
    group1Perm = zeros(size(group1));
    group2Perm = [];
    for k = 1:U
        idx = indexCell{k};
        if isempty(idx)
            continue
        end
        signs = 2 * (rand(rngStream, numel(idx), 1) > 0.5) - 1;
        group1Perm(idx) = group1(idx) .* signs;
    end
end

end
