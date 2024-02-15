function product = celltuples(C)
    %
    % Tuples function of matrices in cell array. Output matrix.
    %
    % @since 1.0.0
    % @param {type} [name] description.
    % @return {type} [name] description.
    % @see dependencies
    %

    D = flip(C);
    Clist = cellfun(@(x) 1:size(x, 1), D, 'UniformOutput', false);
    [D{:}] = ndgrid(Clist{:});
    D = flip(D);
    product = cell2mat(cellfun(@(x, y) x(y, :), C, D, 'UniformOutput', false));
end
