function EXP = expected_condition_counts(Tbefore, excludedSubj, nSubjTotal)
% EXPECTED_CONDITION_COUNTS  Expected number of E/N/R recordings after exclusions.
%   EXP = expected_condition_counts(Tbefore, excludedSubj, nSubjTotal)
%   Full design: nSubjTotal subjects x 2 conditions, i.e. nSubjTotal*2/3 recordings
%   per condition; the recordings of each excluded subject are subtracted.
base = nSubjTotal * 2 / 3;
EXP = struct('E', base, 'N', base, 'R', base);
if isempty(excludedSubj); return; end
s = Tbefore.Subject;
if iscategorical(s) || isstring(s); s = cellstr(s); end
if iscell(s)
    n = cellfun(@(x) str2double(regexprep(char(x), '\D', '')), s);
else
    n = double(s);
end
m = ismember(n, excludedSubj);
rec = unique(Tbefore(m, {'Subject', 'InterventionNumber', 'Condition'}), 'rows');
c = cellstr(rec.Condition);
EXP.E = base - sum(strcmp(c, 'E'));
EXP.N = base - sum(strcmp(c, 'N'));
EXP.R = base - sum(strcmp(c, 'R'));
end
