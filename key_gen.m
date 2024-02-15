function key = key_gen(graph)
    %
    % NLCE unique key generator for different graph.
    %
    % @since 1.0.0
    % @param {type} [name] description.
    % @return {type} [name] description.
    % @see dependencies
    %

    key = logical(adjacency(graph));
    key = num2str(key(:)');
    key = key(key ~= ' ');
    str = ['text2int("', key, '",2)'];
    key = evalin(symengine, str);
    key = sprintf('%s', key);
end
