function T = structrows_to_table(rows)
%STRUCTROWS_TO_TABLE Deterministic scalar-struct-array to table conversion.
%
% This helper exists because struct2table() can interpret cell-valued fields
% differently depending on the exact field shape. Here conversion is explicit:
%   char scalar fields       -> N-by-1 cell array of char
%   numeric/logical scalars  -> N-by-1 numeric/logical array
%   other values             -> N-by-1 cell array
%
% The benchmark row structs are intentionally scalar-valued, so this
% conversion is simple and predictable on MATLAB R2025b and earlier releases.

if isempty(rows)
    T = table();
    return;
end

if ~isstruct(rows)
    error('FiveCodecBenchmark:BadRows','rows must be a struct array.');
end

fields = fieldnames(rows);
n = numel(rows);
T = table();

for k = 1:numel(fields)
    f = fields{k};
    vals = cell(n,1);
    for i = 1:n
        vals{i} = rows(i).(f);
    end

    % Plain char metadata -> cellstr-like table column.
    if all(cellfun(@ischar, vals))
        T.(f) = vals;
        continue;
    end

    % Scalar numeric values -> numeric column.
    isNumScalar = cellfun(@(v) isnumeric(v) && isscalar(v) && isreal(v), vals);
    if all(isNumScalar)
        T.(f) = cellfun(@double, vals);
        continue;
    end

    % Scalar logical values -> logical column.
    isLogScalar = cellfun(@(v) islogical(v) && isscalar(v), vals);
    if all(isLogScalar)
        T.(f) = cellfun(@logical, vals);
        continue;
    end

    % Compatibility fallback: unwrap one scalar cell layer if every value
    % is a 1x1 cell; otherwise retain generic cell values.
    oneCell = cellfun(@(v) iscell(v) && isscalar(v), vals);
    if all(oneCell)
        vals2 = cellfun(@(v) v{1}, vals, 'UniformOutput', false);
        if all(cellfun(@ischar, vals2))
            T.(f) = vals2;
        elseif all(cellfun(@(v) isnumeric(v) && isscalar(v), vals2))
            T.(f) = cellfun(@double, vals2);
        else
            T.(f) = vals2;
        end
    else
        T.(f) = vals;
    end
end
end
