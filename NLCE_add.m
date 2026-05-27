%% Parameters
% Fermi-Hubbard parameters
t = 1; u = [7.9, 13.7, 1.8]; % U12, U13, U23
m = 3;

% Experimental parameters
mfermion = 6*1.66054e-27; %mass of a Li6 atom in kg
trap_laser = 752e-9; % Frequency of the trapping laser (in nm)
h = 6.62607015e-34; % Planck's constant in J-s

% Maximum NLCE order (site-expansion)
order_max = 6;

% Dictionary to map  the pair indices to single integer values.
rev_pair_idx = [1 2; 1 3; 2 3];
%% ED on all graphs up to order order_max.
% Prevents read/write error later.
orderlist=1:order_max;


for order=orderlist
    % Load the graph
    [g, ~] = NLCE_load(order);
    for k = 1:numel(g)
        % Generate unique key for g upto isomorphism.
        obs_list = {"density","doublon","triplon","nearest_pair"};
        item_nums = length(obs_list);
        key = key_gen(g{k});
        l = numnodes(g{k});
        numEdge = numedges(g{k});
        [side1, side2] = findedge(g{k}); %for nearest neighbor correlator

        % Store the matrix elements for observables: <i|O|i> for all eigenstates
        filename = join(['../data/mat_files/ED_mat_files/N=',num2str(m),'/' key],'');
        filename = join([filename," t=", double(t), "u=", double(u), "m=", m, "ED.mat"],' '); 
        
        % Include all the n.n. connections.
        if any(cellfun(@(x) strcmp(x, "nearest_pair"), obs_list))
            obs_list{item_nums + 1} = side1';
            obs_list{item_nums + 2} = side2';
        end

        % Check if the ED data already exists.
        try
            disp(['Try to find ED data for graph ', num2str(key), '...']);
            load(filename);
            varNames = fieldnames(output);
            for i = 1:length(varNames)
                feval(@()assignin('caller', varNames{i}, output.(varNames{i})));
            end
            for idx = 1:numel(obs_list)
                if strcmp(obs_list{idx},"density")
                    if ~exist("nimatrix",'var')
                        delete(filename);
                        load(filename);
                    end
                end
                if strcmp(obs_list{idx},"doublon")
                    if ~exist("domatrix",'var')
                        delete(filename);
                        load(filename);
                    end
                end
                if strcmp(obs_list{idx},"triplon")
                    if ~exist("trimatrix",'var')
                        delete(filename);
                        load(filename);
                    end
                end
                if strcmp(obs_list{idx},"nearest_pair")
                    if ~exist("ninj",'var')
                        delete(filename);
                        load(filename);
                    end
                end
            end
            % Displays the ED run time for the graph g{k}.
            disp(timer_count); 
            disp(['ED data for graph ', num2str(key), ' found.']);
        catch
            disp(['ED data for graph ', num2str(key), ' not found. Start to calculate ED.']);

            % ED function.
            dynamicVars = ED_solver(g{k},t,u,m,filename,obs_list{:});
            % Dynamically assign variables to the workspace
            varNames = fieldnames(dynamicVars);
            for i = 1:length(varNames)
                feval(@()assignin('caller', varNames{i}, dynamicVars.(varNames{i})));
            end
            disp(['ED calculation for graph ', num2str(key), ' finished.']);
            disp(timer_count);
        end
        disp(join([key, 'l=', l, "t=", double(t), "u=", double(u), "m=", m, 'ED measurement start']));
        disp(join(["Measured item :", obs_list{:}]));
    end      
end
    
%% Read the Experimental data
obs_string_arr = {'density','doublon','nearest_pair'}; 
obs_string = 'density';

stringArray{1} = 'density';
for i = 1:numel(u)
    if mod(u(i), 1) == 0  % Check if the number is an integer
        stringArray{end+1} = num2str(u(i));
    else
        stringArray{end+1} = strrep(sprintf('%.1f', u(i)), '.', 'p');
    end
end

% The directory path for the experimental data
directory = join(["../experimental_data_U13_",stringArray{3},"_U12_",stringArray{2},"_U23_",stringArray{4},"/"],"");

% counts and stores the numberof # flavor-combinations for the observable.
col_keys = {'density','doublon','triplon','onsite_pair','nearest_pair'};
col_vals = {m, uint8(m*(m-1)/2), uint8(m*(m-1)*(m-2)/6), uint8(m*(m-1)/2), uint8(m*(m-1)/2)};
col_map = containers.Map(col_keys,col_vals);
col_nums = 0;
for i=1:numel(obs_string_arr)
    col_nums = col_nums + col_map(obs_string_arr{i});
end

% Search for experimental data file and read.
r_file = searchForFile(directory,stringArray);
read_data = readmatrix(r_file,"FileType","text","NumHeaderLines",1);
r = read_data(:,1);
data = zeros(numel(r),col_nums);
errors = zeros(numel(r),col_nums);
data(:,1:m) = read_data(:,2:m+1);
errors(:,1:m) = read_data(:,m+2:m+4);
column_number = m; % counter to assign the correct column;

% Normalize data with the highest value.
norm_data = data./max(abs(data));

% Extract HTE0 fit from file (mu0's, T, lattice confinement and tunneling)
[dirPath,~,~] = fileparts(r_file);
parameter_file_path = fullfile(dirPath,'Parameters.txt');
[mu0_exp, T_exp, omega, tunneling, N] = find_exp_params(parameter_file_path,u,m);

%% Single observable atomic limit fit
obs_string = 'density'; % Any of the other available observables are acceptable.
indices = find(r<42);
xData = r(indices);

% Known system parameters
otherParams = {t,u,m,1,tunneling, omega, {'density'}, mfermion, trap_laser, h};
yData = data(indices,1:m);

% atomic limit fit. statistics to makes sure, not getting stuck in local
% minimum.
% Fitting parameters: Temperature, chemical potential at the ceneter of the
% trap for the three flavors.
TGuess = rand(1); %temperature
muguess = rand(3,1); % center-of-trap chemical potentials.
fit_atomic_limit = zeros(4);

initialGuess = [TGuess, muguess(:)'];  
options = optimoptions('lsqcurvefit' ,'OptimalityTolerance', 1e-6,'MaxIterations',10000,'FunctionTolerance',1e-6, 'Algorithm','trust-region-reflective','Display','off');
[fittingParams, resnorm,residual,exitflag,output,lambda,J] = lsqcurvefit(@(fitParams,x) obs_vs_r_fitting_function(x,fitParams,otherParams), initialGuess, xData, yData,[1,-10,-10,-10],[5,10,10,10],options);
fit_atomic_limit = fittingParams;
covariance = inv(J' * J)*resnorm/(numel(yData) - numel(fittingParams));
std_err = full(sqrt(diag(covariance)));


ustr = regexprep(num2str(u),'\s+', ', ');
ustr = strrep(ustr, '.', 'p');
fitparam_file = join(["../data/csv_files/N=3/", "u=", ustr  "_fit_vals.txt"],'');
fileID = fopen(fitparam_file,'w');
fprintf(fileID,"Atomic Limit Fit\n");
fprintf(fileID, "T/t = %.4f +/- %.4f\nmu0/t = ", fit_atomic_limit(1),std_err(1));
for i=2:4
    fprintf(fileID, "%.4f +/- %.4f\t", fit_atomic_limit(i), std_err(i));
end
fprintf(fileID,repmat('\n',1,4));
fclose(fileID);

%% Bootstrapping the fit to the desired NLCE order
obs_string_arr = {'density'};
yData = data(indices,1:m);
fitGuess = fit_atomic_limit;
for fit_order = 2:order_max
    otherParams = {t,u,m,fit_order,tunneling, omega, obs_string_arr, mfermion, trap_laser, h};
    [fittingParams, resnorm, residual, exitflag, output, lambda, J] = lsqcurvefit(@(fitParams,x) obs_vs_r_fitting_function(x,fitParams,otherParams), fitGuess, xData, yData,[.5,-10,-10,-10],[5,10,10,10],options);
    fitGuess = fittingParams;
end
covariance = inv(J' * J)*resnorm/(numel(yData) - numel(fittingParams));
std_err = full(sqrt(diag(covariance)));
fitparam_file = join(["../data/csv_files/N=3/", "u=", ustr  "_fit_vals.txt"],'');
fileID = fopen(fitparam_file,'a'); % This will append to the original file! Use this in conjunction with the atomic limit fit.
fprintf(fileID,"NLCE order: %d\n",order_max);
fprintf(fileID, "T/t = %.4f +/- %.4f\nmu0/t = ", fitGuess(1), std_err(1));
for i=2:4
    fprintf(fileID, "%.4f +/- %.4f\t", fitGuess(i), std_err(i));
end
fclose(fileID);

%% mu-T grid of from the fitted values.
muq = zeros(numel(r),m);
for i=1:m
    muq(:,i) =  fittingParams(i+1) - .5*(mfermion)*(2*pi*omega)^2*(r*trap_laser).^2/(h*tunneling); 
    % 752 nm is the lattice spacing distance.
end

Tarray = [fittingParams(1)];
mu_file = join(['../data/csv_files/N=',num2str(m), '/u=', ustr, '_mus.csv'],'');
T_file = join(['../data/csv_files/N=',num2str(m),'/u=', ustr, '_T.csv'],'');
writematrix(muq,mu_file);
writematrix(Tarray,T_file);
clear("mu_file","T_file");

%% Generate the NLCE sum for all observables up to order_max.

orderlist = 1:order_max;

% Dictionary of observables to be computed
obs_list = {"density","doublon","triplon","onsite_pair","nearest_pair"};

NLCE_sum(t,u,m,obs_list,orderlist,muq,Tarray,true); % 4D double [m,mu,T,order]
% Keep the last entry to be true for saving data.

%% Functions
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
                            densityaccum(:,spin_idx) = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, rho_matrix);
                            % UNCOMMENT FOR SAME SPIN CORRELATIONS
                            % if l>1
                            %     ninj_matrix = cellfun(@(x) x{spin_idx+uint8(m*(m-1)/2)}, ninj,'UniformOutput',false);
                            %     nn_correlator_accum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, ninj_matrix);
                            %     nearest_pair_accum(uint8(m*(m-1)/2)+spin_idx) = sum(nn_correlator_accum(1:numEdge)- densityaccum(side1,spin_idx).*densityaccum(side2,spin_idx));
                            % end
                        end
                        doaccum = zeros(1,uint8(m*(m-1)/2));
                        for pair_idx = 1:uint8(m*(m-1)/2)
                            pair_matrix = cellfun(@(x) x{pair_idx},domatrix,'UniformOutput',false);
                            doaccum(pair_idx) = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, pair_matrix);
                            if l>1
                                ninj_matrix = cellfun(@(x) x{pair_idx}, ninj,'UniformOutput',false);
                                nn_correlator_accum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, ninj_matrix);
                                nearest_pair_accum(pair_idx) = sum(nn_correlator_accum(1:numEdge)- densityaccum(side1,rev_pair_idx(pair_idx,1)).*densityaccum(side2,rev_pair_idx(pair_idx,2)));
                            end
                        end
                        
                        for triplon_idx = 1:uint8(m*(m-1)*(m-2)/6)
                            triplon_matrix = cellfun(@(x) x{triplon_idx},trimatrix,'UniformOutput',false);
                            triaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, triplon_matrix);
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
                            densityaccum(:,spin_idx) = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, rho_matrix);
                        end
                        temp_list(i,column_num+1:column_num + m) = sum(densityaccum,1);
                        column_num = column_num + m;
                    elseif strcmp(obs_string, 'doublon')
                        doaccum = zeros(1,uint8(m*(m-1)/2));
                        for pair_idx = 1:uint8(m*(m-1)/2)
                            pair_matrix = cellfun(@(x) x{pair_idx},domatrix,'UniformOutput',false);
                            doaccum(pair_idx) = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, pair_matrix);
                        end
                        temp_list(i,column_num + 1: column_num + uint8(m*(m-1)/2)) = doaccum(:)';
                        column_num = column_num + uint8(m*(m-1)/2);
                    elseif strcmp(obs_string, 'triplon')
                        triaccum = zeros(1,uint8(m*(m-1)*(m-2)/6));
                        for tri_idx = 1:uint8(m*(m-1)*(m-2)/6)
                            triplon_matrix = cellfun(@(x) x{tri_idx}, trimatrix, 'UniformOutput',false);
                            triaccum = LDA_measure(spectra,deltamu, T, nsigmapermute, testn, triplon_matrix);
                        end
                        temp_list(i,column_num + 1:column_num + uint8(m*(m-1)*(m-2)/6)) = triaccum;
                        column_num = column_num + uint8(m*(m-1)*(m-2)/6);
                    elseif strcmp(obs_string, 'nearest_pair')
                        nearest_pair_accum = zeros(1,uint8(m*(m-1)/2));
                        for pair_idx = 1:uint8(m*(m-1)/2)
                            if l>1
                                ninj_matrix = cellfun(@(x) x{pair_idx}, ninj,'UniformOutput',false);
                                nn_correlator_accum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, ninj_matrix);
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

function [mu0, T, omega, tunneling, N] = find_exp_params(filepath, u, m)
    
    % Open the file for reading
    fileID = fopen(filepath, 'r');
    
    % Initialize variables to store the result
    mu0 = []; T = NaN; omega = NaN; tunneling = NaN; N = [];

    % Define the string pattern after which "mu1" should occur
    exp_u = [u(2) u(1) u(3)]; % experiment data has U13 U12 U23.
    pattern = sprintf("%.1f ",exp_u);
    pattern = strtrim(pattern);
   % Read each line of the file
    found_pattern = false;
    tline = fgetl(fileID);
    lineCounter = 0;
    while ~found_pattern && ischar(tline)
        
        % Search for the specific pattern in the line
        idx_pattern = strfind(tline, pattern);
        if ~isempty(idx_pattern)
            found_pattern = true; % Set flag to true if the pattern is found
        end
        % Read the next line
        tline = fgetl(fileID);
        lineCounter = lineCounter + 1;
    end
    
    % Once the pattern is found, search for the next occurrence of "mu_i"
    % and "T"
    for i=1:m+1
        while ischar(tline) && found_pattern
            % Search for "mu_i" and "T" in the line
            idx_T = strfind(tline, "T:");
            idx_mu = strfind(tline, sprintf('mu%d:',i-1));          
            if i==1
                if ~isempty(idx_T)
                    T_str = tline(idx_T(1) + 3:end);
                    T = str2double(T_str);
                    break;
                end
            end
            if ~isempty(idx_mu)
                % If "mu1" is found, extract the float value
                mu_str = tline(idx_mu(1) + 5:end); % Skip "mu1: "
                mu0(end+1) = str2double(mu_str); % Convert to float
                break; % Exit the loop once found
            end
            % Read the next line
            tline = fgetl(fileID);
        end
    end

    %  Look for "t=<number>"
    % Calculate the line number nLinesBefore the found string
    targetLineNumber =  lineCounter - 3;

    % Go to the target line
    fseek(fileID, 0, 'bof');
    for i = 1:targetLineNumber-1
        fgetl(fileID);
    end
    
    prevLine = fgetl(fileID);
    % Look for the number in the previous line
    if contains(prevLine, 't=')
        tunneling = str2double(regexp(prevLine, 't=(\d+\.?\d*)', 'tokens', 'once'));
    end
    
    %Look for "x <number>"
    nextLine = fgetl(fileID);
    if contains(nextLine,"x ")
        omega = str2double(regexp(nextLine, 'x (\d+\.?\d*)', 'tokens', 'once'));
    end

    %Look for "(n1,n2,n3) = <number> "
    nextLine = fgetl(fileID);
    if contains(nextLine,"(n1,n2,n3)")
        % Define a regular expression pattern to extract the numbers
        pattern = '\(([^)]+)\)'; % Pattern to match content inside parentheses
        
        % Extract the content inside parentheses
        matches = regexp(nextLine, pattern, 'match');
        
        % Extract the numbers from the second match (right-hand side of the equation)
        numbers_str = matches{2};  % Get the second matched group, which contains the numbers
        
        % Remove the parentheses
        numbers_str = numbers_str(2:end-1);
        
        % Split the string by commas to get individual number strings
        number_cells = strsplit(numbers_str, ',');
        
        % Convert the cell array of strings to an array of numbers
        N = str2double(number_cells);
    end
    
    % Close the file
    fclose(fileID);
end


function foundFile = searchForFile(directory, searchStrings)
    % Initialize the output
    foundFile = '';

    % Get a list of all files in the current directory
    files = dir(directory);

    % Loop through each file in the current directory
    for j = 1:length(files)
        % Get the filename
        filename = files(j).name;

        % Check if the filename contains all the search strings
        containsAllStrings = all(cellfun(@(str) contains(filename, str), searchStrings));

        % If a match is found, set the output and return
        if containsAllStrings
            foundFile = fullfile(directory, filename);
            return;
        end
    end
    
    % If no match is found, display a message
    if isempty(foundFile)
        disp('No file found matching the given search strings.');
    end
end