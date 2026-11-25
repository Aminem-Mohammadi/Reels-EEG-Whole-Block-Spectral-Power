function [T, nRowsExcluded, subjectsExcluded] = apply_run_mode_filter(T, RUN_MODE, EXCLUDE_LIST)
% APPLY_RUN_MODE_FILTER  Remove the cleaned-sample exclusions from a table or struct array.
%   [T, nRowsExcluded, subjectsExcluded] = apply_run_mode_filter(T, RUN_MODE, EXCLUDE_LIST)
%   Does nothing unless RUN_MODE is 'sensitivity_excl_noisy'. T is a table with a
%   'Subject' column or a struct array with a 'Subject' field; subject IDs may be
%   numeric or text (digits are extracted, e.g. 'sub07' -> 7).

nRowsExcluded = 0;
subjectsExcluded = [];
if ~strcmp(RUN_MODE, 'sensitivity_excl_noisy') || isempty(EXCLUDE_LIST) || isempty(T)
    return;
end
if isstruct(T)   % struct array (e.g. psd_data entries) with a Subject field
    if ~isfield(T, 'Subject')
        error('RunMode:NoSubjectField', 'apply_run_mode_filter: struct array has no Subject field.');
    end
    n = subject_to_number({T.Subject});
    mask = ismember(n, EXCLUDE_LIST);
    subjectsExcluded = unique(n(mask))';
    nRowsExcluded = sum(mask);
    T = T(~mask);
    return;
end
if ~ismember('Subject', T.Properties.VariableNames)
    error('RunMode:NoSubjectColumn', 'apply_run_mode_filter: table has no Subject column.');
end
n = subject_to_number(T.Subject);
mask = ismember(n, EXCLUDE_LIST);
subjectsExcluded = unique(n(mask))';
nRowsExcluded = sum(mask);
T = T(~mask, :);
end

function n = subject_to_number(s)
if iscategorical(s) || isstring(s)
    s = cellstr(s);
end
if isnumeric(s)
    n = double(s(:));
elseif iscell(s)
    n = cellfun(@(x) local_one(x), s(:));
else
    n = local_one(s);
end
end

function v = local_one(x)
if isnumeric(x)
    v = double(x);
else
    v = str2double(regexprep(char(x), '\D', ''));
end
end
