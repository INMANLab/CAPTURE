function s = cohortMeanStatistic(participantValues)
%COHORTMEANSTATISTIC Cohort-level summary: mean across participants.
%
%   Formula
%   -------
%   Let m_p be the scalar summary for participant p (e.g. mean of trials
%   within that participant). With P participants,
%
%       s = (1 / P_valid) * sum_p m_p
%
%   where only finite (non-NaN) m_p are included (P_valid = count of finite
%   entries). This matches MATLAB mean(..., 'omitnan') over the P-vector.
%
%   If all entries are NaN, s is NaN.

    arguments
        participantValues (:,1) double
    end

    s = mean(participantValues, 'omitnan');
end
