function obs_arr = obs_vs_r_fitting_function(x, fitParams, otherParams)
% Returns normalized arrays of the observable
   
    T = fitParams(1);
    mu0 = fitParams(2:4);
    rev_pair_idx = [1,2;1,3;2,3];

    t = otherParams{1}; u = otherParams{2}; m = otherParams{3};
    order_max = otherParams{4}; tunneling = otherParams{5}; omega = otherParams{6};
    obs_string_arr = otherParams{7};
    mfermion = otherParams{8}; trap_laser = otherParams{9}; h = otherParams{10};
    

    column_count = 0;
    col_keys = {'density','doublon','triplon','nearest_pair'};
    col_vals = {m,uint8(m*(m-1)/2),uint8(m*(m-1)*(m-2)/6),uint8(m*(m-1)/2)};
    col_map = containers.Map(col_keys,col_vals);
    for i=1:numel(obs_string_arr)
        obs_string = obs_string_arr{i};
        column_count = column_count + col_map(obs_string);
    end
    
    mu = zeros(numel(x),m);
    for i = 1:m
        mu(:,i) = mu0(i) - .5*(mfermion)*(2*pi*omega)^2*(x*trap_laser).^2/(h*tunneling); 
    end
    obs_arr = zeros();
    for order=1:order_max
        [g, coefficient] = NLCE_load(order);
        for k = 1:numel(g)
            key = key_gen(g{k});
            l = numnodes(g{k});
            filename = join(['../data/mat_files/ED_mat_files/N=',num2str(m),'/', key],'');
            filename = join([filename," t=", double(t), "u=", double(u), "m=", m, "ED.mat"],' '); % name of file that saves the output
            if ismember('nearest_pair',{obs_string_arr{1:numel(otherParams{9})}})
                [side1,side2] = findedge(g{k});
                obs_string_arr{column_count + 1} = side1';
                obs_string_arr{column_count + 2} = side2';
                numEdge = numedges(g{k});
            end
            disp(['Try to find ED data for graph ', num2str(key), '...']);
            % Try to load ED result of graph. Run ED if not found.
            try
                load(filename); 
                varNames = fieldnames(output);
                for i = 1:length(varNames)
                    feval(@()assignin('caller', varNames{i}, output.(varNames{i})));
                end
                disp(['ED data for graph ',num2str(key), ' found.']);
            catch
                disp(['ED data for graph ', num2str(key), ' not found. Start to calculate ED.']);
                dynamicVars = ED_solver(g{k},t,u,m,filename,obs_string_arr{:});
                % Dynamically assign variables to the workspace
                varNames = fieldnames(dynamicVars);
                for i = 1:length(varNames)
                    feval(@()assignin('caller', varNames{i}, dynamicVars.(varNames{i})));
                end
                disp(['ED calculation for graph ', num2str(key), ' finished.']);
                disp(timer_count);
            end
            
            temp_list = zeros(numel(x),column_count);
            % For every temperature, perform NLCE sum
            for i = 1:numel(x)
                % Fix mu=0 to be half-filling
                deltamu =  - mu(i,:); 
                column_num = 0;
                for j=1:numel(obs_string_arr)
                    obs_string = obs_string_arr{j};
                    if strcmp(obs_string, 'density')
                        densityaccum = zeros(l,m);
                        for spin_idx = 1:m
                            rho_matrix = cellfun(@(x) x{spin_idx},nimatrix,'UniformOutput',false);
                            densityaccum(:,spin_idx) = thermal_average(spectra, deltamu, T, nsigmapermute, testn, rho_matrix);
                        end
                        temp_list(i,column_num+1:column_num + m) = sum(densityaccum,1);
                        column_num = column_num + m;
                    elseif strcmp(obs_string, 'doublon')
                        doaccum = zeros(1,uint8(m*(m-1)/2));
                        for pair_idx = 1:uint8(m*(m-1)/2)
                            pair_matrix = cellfun(@(x) x{pair_idx},domatrix,'UniformOutput',false);
                            doaccum(pair_idx) = thermal_average(spectra, deltamu, T, nsigmapermute, testn, pair_matrix);
                        end
                        temp_list(i,column_num + 1: column_num + uint8(m*(m-1)/2)) = doaccum(:)';
                        column_num = column_num + uint8(m*(m-1)/2);
                    elseif strcmp(obs_string, 'triplon')
                        triaccum = zeros(1,uint8(m*(m-1)*(m-2)/6));
                        for tri_idx = 1:uint8(m*(m-1)*(m-2)/6)
                            triplon_matrix = cellfun(@(x) x{tri_idx}, trimatrix, 'UniformOutput',false);
                            triaccum = thermal_average(spectra,deltamu, T, nsigmapermute, testn, triplon_matrix);
                        end
                        temp_list(i,column_num + 1:column_num + uint8(m*(m-1)*(m-2)/6)) = triaccum;
                        column_num = column_num + uint8(m*(m-1)*(m-2)/6);
                    elseif strcmp(obs_string, 'nearest_pair')
                        nearest_pair_accum = zeros(1,uint8(m*(m-1)/2));
                        for pair_idx = 1:uint8(m*(m-1)/2)
                            if l>1
                                ninj_matrix = cellfun(@(x) x{pair_idx}, ninj,'UniformOutput',false);
                                nn_correlator_accum = thermal_average(spectra, deltamu, T, nsigmapermute, testn, ninj_matrix);
                                nearest_pair_accum(pair_idx) = sum(nn_correlator_accum(1:numEdge)- densityaccum(side1,rev_pair_idx(pair_idx,1)).*densityaccum(side2,rev_pair_idx(pair_idx,2)));
                            end
                        end
                        temp_list(i, column_num + 1 : column_num + uint8(m*(m-1)/2) ) = nearest_pair_accum(:)';
                    end
                end
            end
            obs_arr = obs_arr + coefficient(k)*temp_list; %nlce sum
            disp(join([key, 'l=', l, "t=", double(t), "u=", double(u), "m=", m, 'ED measurement finish']));
        end
    end
end
