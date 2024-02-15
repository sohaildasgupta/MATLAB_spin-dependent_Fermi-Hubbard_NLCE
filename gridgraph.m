function g = gridgraph(varargin)
    %
    % g = gridgraph(sz1, sz2, sz3, option)
    %
    % Only support up to 3-D grid.
    % option specified Boudary Codtion : 'obc' / 'pbc'
    %
    % @since 1.0.0
    % @param {type} [name] description.
    % @return {type} [name] description.
    % @see dependencies
    %

    input = [varargin{1:end-1}];
    option = varargin{end};
    dim = numel(input);
    l = prod(input);

    if l == 1
        g = graph(1);
    else

        if dim - 1
            nodeidx = reshape(1:l, input);
        else
            nodeidx = (1:l)';
        end

        s = []; t = [];

        if strcmp(option, 'pbc')

            for idx = 1:dim
                bond = circshift(nodeidx, 1, idx);
                s = [s; nodeidx(:)];
                t = [t; bond(:)];
            end

        else

            for idx = 1:dim
                dimorder = circshift(1:dim, idx - 1);

                if dim - 1
                    permnodeidx = permute(nodeidx, dimorder);
                else
                    permnodeidx = nodeidx;
                end

                startnode = permnodeidx(1:end - 1, :, :);
                endnode = permnodeidx(2:end, :, :);
                s = [s; startnode(:)];
                t = [t; endnode(:)];
            end

        end

        g = graph(s, t);
    end

end
