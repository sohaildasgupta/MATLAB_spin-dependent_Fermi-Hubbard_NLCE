function [g, coefficient] = NLCE_load(n)
    %
    % [graph, coefficient] = NLCE_load(n)
    %
    % @since 1.0.0
    % @param {type} [name] description.
    % @return {type} [name] description.
    % @see dependencies
    %

    location = 'NLCE/';

    connection = importdata([location, 'NLCE 2D graphs/graphsSimplified', num2str(n), '.txt']);
    connection = cellfun(@(y) cell2mat(cellfun(@(x) cell2mat(x)', y, 'un', 0))', cellfun(@(x) eval(x), connection, 'un', 0), 'un', 0);
    g = {};
    if numel(connection{1})==0
        g{1} = graph(1);
        g = [g, cellfun(@(x) graph(x(:, 1), x(:, 2)), {connection{2:end}}, 'un', 0)];
    else
        g = [g, cellfun(@(x) graph(x(:, 1), x(:, 2)), {connection{1:end}}, 'un', 0)];
    end
    
    coefficient = importdata([location, 'NLCE 2D coefficients/coefficientsOfGraphs', num2str(n), '.txt']);
end  
