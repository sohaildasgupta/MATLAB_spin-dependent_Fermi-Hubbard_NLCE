function varOutputs = NLCE_sum(t,u,m,observables,orderlist,muq,Tarray,save_files)
% Compute the NLCE sum for observables for given orders, chemical potentials and temperatures.
% INPUTS
%   t: tunneling 
%   u: Hubbard interaction. [u12, u13, u23]
%   m: # spin flavors
%   observables: List of observables to compute. 
%               Options: "density", "doublon", "triplon", "nearest_pair"
%   orderlist: NLCE orders.
%   muq: Chemical potentials 
%   Tarray: Temperatures
%   save_files: Boolean- True to save to files.
%
%   OUTPUTS
%   varoutputs: NLCE_<observable> - 3D array [# mus, # Ts, # orders]. 
%               NLCE sum of the thermal averages of the observable.
    
    rev_pair_idx = [1,2;1,3;2,3];
    itemlist = cellfun(@(x) x, observables, 'UniformOutput',false);
    item_nums = length(itemlist);
    for i=1:numel(itemlist)
        if strcmp(itemlist{i},"density")
            NLCE_density = zeros(m,size(muq,1),numel(Tarray),numel(orderlist));
        end
        if strcmp(itemlist{i},"doublon")
            NLCE_doublon = zeros(uint8(m*(m-1)/2),size(muq,1),numel(Tarray),numel(orderlist));
        end
        if strcmp(itemlist{i},"triplon")
            NLCE_triplon = zeros(size(muq,1),numel(Tarray),numel(orderlist));
        end
        if strcmp(itemlist{i},"onsite_pair")
            NLCE_onsite_pair = zeros(uint8(m*(m-1)/2),size(muq,1),numel(Tarray),numel(orderlist));
        end
        if strcmp(itemlist{i},"nearest_pair")
            NLCE_nearest_pair = zeros(m+uint8(m*(m-1)/2),size(muq,1),numel(Tarray),numel(orderlist));
        end
    end
    CTlist = zeros();
    for order=orderlist
        % order=2; % Just for check
        [g, coefficient] = NLCE_load(order);
        for k = 1:numel(g)
            % k=2; % for checking
            key = key_gen(g{k});
            l = numnodes(g{k});
            numEdge = numedges(g{k});
            [side1, side2] = findedge(g{k}); %for nearest neighbor correlator
            % side1 = [side1' 1:l]; side2 = [side2' 1:l]; %on-site correlator
            filename = join(['../data/mat_files/ED_mat_files/N=',num2str(m),'/' key],'');
            filename = join([filename," t=", double(t), "u=", double(u), "m=", m, "ED.mat"],' '); % name of file that saves the output
            if any(cellfun(@(x) strcmp(x, "nearest_pair"), itemlist))
                itemlist{item_nums + 1} = side1';
                itemlist{item_nums + 2} = side2';
            end
            try
                disp(['Try to find ED data for graph ', num2str(key), '...']);
                % clear('output');
                load(filename);
                varNames = fieldnames(output);
                for i = 1:length(varNames)
                    feval(@()assignin('caller', varNames{i}, output.(varNames{i})));
                end
                for idx = 1:numel(itemlist)
                    if strcmp(itemlist{idx},"density")
                        if ~exist("nimatrix",'var')
                            delete(filename);
                            load(filename);
                        end
                    end
                    if strcmp(itemlist{idx},"doublon")
                        if ~exist("domatrix",'var')
                            delete(filename);
                            load(filename);
                        end
                    end
                    if strcmp(itemlist{idx},"triplon")
                        if ~exist("trimatrix",'var')
                            delete(filename);
                            load(filename);
                        end
                    end
                    if strcmp(itemlist{idx},"nearest_pair")
                        if ~exist("ninj",'var')
                            delete(filename);
                            load(filename);
                        end
                    end
                end
                disp(timer_count);
                disp(['ED data for graph ', num2str(key), ' found.']);
            catch
                disp(['ED data for graph ', num2str(key), ' not found. Start to calculate ED.']);
    
                dynamicVars = ED_solver(g{k},t,u,m,filename,itemlist{:});
                % Dynamically assign variables to the workspace
                varNames = fieldnames(dynamicVars);
                for i = 1:length(varNames)
                    feval(@()assignin('caller', varNames{i}, dynamicVars.(varNames{i})));
                end
                disp(['ED calculation for graph ', num2str(key), ' finished.']);
                disp(timer_count);
            end
    
            disp(join([key, 'l=', l, "t=", double(t), "u=", double(u), "m=", m, 'ED measurement start']));
            disp(join(["Measured item :", itemlist{:}]));
    
            try
                nlce_path = join(["../data/mat_files/NLCE_mat_files/N=",num2str(m),"/"],'');
                nlce_filename = join([key," t=", double(t), "u=", double(u), "m=", m, "LDA"],' ');
                nlce_filename = strrep(nlce_filename,'.','p');
                nlce = join([nlce_path,nlce_filename,".mat"],"");
                load(nlce);
                
                %WARNING : This deletes previous data. Need a better check.
                if ~isequal(mulist, muq)
                    disp("mu-list do not match, deleting nlce graph data...");
                    delete(nlce);
                    load(nlce);
                end
                if ~isequal(Tlist, Tarray)
                    disp("T-list do not match, deleting nlce graph data...");
                    delete(nlce);
                    load(nlce);
                end
            catch
                alist = [];
                for T = Tarray
                    for mu_idx = 1:size(muq,1)
                        mu0 = muq(mu_idx,:);
                        deltamu =  - mu0; % REWRITE To include different mus
                        densityaccum = zeros(l,m);
                        nearest_pair_accum = zeros(1,m+uint8(m*(m-1)/2));
                        for spin_idx = 1:m
                            rho_matrix = cellfun(@(x) x{spin_idx},nimatrix,'UniformOutput',false);
                            densityaccum(:,spin_idx) = thermal_average(spectra, deltamu, T, nsigmapermute, testn, rho_matrix);
                            % UNCOMMENT FOR SAME SPIN CORRELATIONS
                            % if l>1
                            %     ninj_matrix = cellfun(@(x) x{spin_idx+uint8(m*(m-1)/2)}, ninj,'UniformOutput',false);
                            %     nn_correlator_accum = thermal_average(spectra, deltamu, T, nsigmapermute, testn, ninj_matrix);
                            %     nearest_pair_accum(uint8(m*(m-1)/2)+spin_idx) = sum(nn_correlator_accum(1:numEdge)- densityaccum(side1,spin_idx).*densityaccum(side2,spin_idx));
                            % end
                        end
                        doaccum = zeros(1,uint8(m*(m-1)/2));
                        for pair_idx = 1:uint8(m*(m-1)/2)
                            pair_matrix = cellfun(@(x) x{pair_idx},domatrix,'UniformOutput',false);
                            doaccum(pair_idx) = thermal_average(spectra, deltamu, T, nsigmapermute, testn, pair_matrix);
                            if l>1
                                ninj_matrix = cellfun(@(x) x{pair_idx}, ninj,'UniformOutput',false);
                                nn_correlator_accum = thermal_average(spectra, deltamu, T, nsigmapermute, testn, ninj_matrix);
                                nearest_pair_accum(pair_idx) = sum(nn_correlator_accum(1:numEdge)- densityaccum(side1,rev_pair_idx(pair_idx,1)).*densityaccum(side2,rev_pair_idx(pair_idx,2)));
                            end
                        end
                        
                        for triplon_idx = 1:uint8(m*(m-1)*(m-2)/6)
                            triplon_matrix = cellfun(@(x) x{triplon_idx},trimatrix,'UniformOutput',false);
                            triaccum = thermal_average(spectra, deltamu, T, nsigmapermute, testn, triplon_matrix);
                        end
                        alist(end + 1, :) = [sum(densityaccum,1) doaccum(:)' triaccum nearest_pair_accum(:)'];
                    end
                end
                nlce_path = join(["../data/mat_files/NLCE_mat_files/N=",num2str(m),"/"],'');
                nlce_filename = join([key," t=", double(t), "u=", double(u), "m=", m, "LDA"],' ');
                nlce_filename = strrep(nlce_filename,'.','p');
                nlce = join([nlce_path,nlce_filename,".mat"],"");
                mulist = muq; %To ensure the correct mus are used
                Tlist = Tarray;
                save(nlce,'alist','mulist','Tlist');
            end
            CTlist = CTlist + coefficient(k)*alist; %nlce sum
            disp(join([key, 'l=', l, "t=", double(t), "u=", double(u), "m=", m, 'ED measurement finish']));
        end

        % Convert the u array to string
        if save_files 
            ustr = regexprep(num2str(u),'\s+', ', ');
            ustr = strrep(ustr, '.', 'p');
        end

        % Return the desired data.
        if exist("NLCE_density",'var') == 1
            for idx = 1:m
                NLCE_density(idx,:,:,order) = reshape(CTlist(:, idx),[size(muq,1),numel(Tarray)])';
                if save_files
                    writematrix(NLCE_density(idx,:,:,order)', join(['../data/csv_files/N=' num2str(m) '/NLCE_order=' num2str(order) '_u=' ustr  '_density_flavor=' num2str(idx) '.csv'],''));
                end
            end
        end

        if exist("NLCE_doublon",'var') == 1
            for idx = 1:uint8(m*(m-1)/2)
                NLCE_doublon(idx,:,:,order) = reshape(CTlist(:, m+idx),[size(muq,1),numel(Tarray)])';
                if save_files
                    writematrix(NLCE_doublon(idx,:,:,order)', join(['../data/csv_files/N=' num2str(m) '/NLCE_order=' num2str(order) '_u=' ustr  '_doublon_pair=' strrep(num2str(rev_pair_idx(idx,:)),' ','') '.csv'],''));
                end
            end
        end

        if exist("NLCE_triplon",'var') == 1
            NLCE_triplon(:,:,order) = reshape(CTlist(:, m+uint8(m*(m-1)/2)+1),[size(muq,1),numel(Tarray)])';
            if save_files
                writematrix(NLCE_triplon(:,:,order), join(['../data/csv_files/N=' num2str(m) '/NLCE_order=' num2str(order) '_u=' ustr '_triplon.csv'],''));
            end
        end
        
        if exist("NLCE_onsite_pair",'var') == 1
            for idx = 1:uint8(m*(m-1)/2)
                NLCE_onsite_pair(idx,:,:,order) = NLCE_doublon(idx,:,:,order)-  NLCE_density(rev_pair_idx(idx,1),:,:,order).*NLCE_density(rev_pair_idx(idx,2),:,:,order);
                if save_files
                    writematrix(NLCE_onsite_pair(idx,:,:,order)', join(['../data/csv_files/N=' num2str(m) '/NLCE_order=' num2str(order) '_u=' ustr  '_onsite_pair=' strrep(num2str(rev_pair_idx(idx,:)),' ','') '.csv'],''));
                end
            end
        end

        if exist("NLCE_nearest_pair",'var') == 1
            if l>1
                for idx = 1:m+uint8(m*(m-1)/2)
                    NLCE_nearest_pair(idx,:,:,order) = reshape(CTlist(:, m+uint8(m*(m-1)/2)+1+idx),[size(muq,1),numel(Tarray)])';% - NLCE_density(rev_pair_idx(idx,1),:,:,order).*NLCE_density(rev_pair_idx(idx,2),:,:,order);
                    if save_files
                        if idx<=m
                            pair_string = strrep(num2str(rev_pair_idx(idx,:)),' ','');
                        else
                            pair_string = sprintf("%i%i",[idx-m,idx-m]);
                        end
                        writematrix(NLCE_nearest_pair(idx,:,:,order)', join(['../data/csv_files/N=' num2str(m) '/NLCE_order=' num2str(order) '_u=' ustr  '_nearest_pair=' pair_string '.csv'],''));
                    end
                end
            else
                NLCE_nearest_pair(:,:,:,order) = 0;
            end
            
        end
    end
    varOutputs = struct();
    if exist("NLCE_density",'var')==1
        varOutputs.density = NLCE_density;
    end
    if exist("NLCE_doublon",'var')==1
        varOutputs.doublon = NLCE_doublon;
    end
    if exist("NLCE_triplon",'var')==1
        varOutputs.triplon = NLCE_triplon;
    end
    if exist("NLCE_onsite_pair",'var')==1
        varOutputs.onsite_pair = NLCE_onsite_pair;
    end
    if exist("NLCE_nearest_pair",'var')==1
        varOutputs.nearest_pair = NLCE_nearest_pair;
    end
end
