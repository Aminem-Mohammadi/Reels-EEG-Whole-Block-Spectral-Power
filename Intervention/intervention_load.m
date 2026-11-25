function S = intervention_load(varargin)
% INTERVENTION_LOAD  load() plus two compatibility columns.
%   S = intervention_load(file, var1, var2, ...) calls load() and then adds, to every
%   loaded table or struct array that has an InterventionNumber column/field,
%   TaskType = 'intervention' and BlockNumber = InterventionNumber (column names
%   expected by the shared analysis code).

S = load(varargin{:});
fn = fieldnames(S);
for k = 1:numel(fn)
    v = S.(fn{k});
    if istable(v) && ismember('InterventionNumber', v.Properties.VariableNames)
        if ~ismember('TaskType', v.Properties.VariableNames)
            v.TaskType = repmat({'intervention'}, height(v), 1);
        end
        if ~ismember('BlockNumber', v.Properties.VariableNames)
            v.BlockNumber = double(v.InterventionNumber);
        end
        S.(fn{k}) = v;
    elseif isstruct(v) && ~isempty(v) && isfield(v, 'InterventionNumber')
        for i = 1:numel(v)
            if ~isfield(v, 'TaskType') || isempty(v(i).TaskType)
                v(i).TaskType = 'intervention';
            end
            if ~isfield(v, 'BlockNumber') || isempty(v(i).BlockNumber)
                v(i).BlockNumber = double(v(i).InterventionNumber);
            end
        end
        S.(fn{k}) = v;
    end
end
end
