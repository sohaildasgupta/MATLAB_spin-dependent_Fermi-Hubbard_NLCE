% Performs the NLCE sum.
% Calls NLCE_load to read the graphs and their corresponding coefficients.
% Calls ED_NLCE to perform the ED on each graph.
% Adds the properties together multiplying by the weight of each graph

%% ED Parameters 
feature('numcores')
t = 1; u = [32.6 87.0 15.5];
m = 3;
n = -1;
dno = -1;

%Dictionary to map  the pair indices to single integer values.
rev_pair_idx = [1 2; 1 3; 2 3];

%% Obtaining the experimental data
file_path = "../experimental_data/";
% Step 1: Convert float array to string array with "." replaced by "p"
stringArray = arrayfun(@(x) strrep(sprintf('%.1f', x), '.', 'p'), u, 'UniformOutput', false);

% Step 2: Loop through directory to check for filenames
directory = '../experimental_data'; % Specify your directory path
fileNames = dir(directory);
for i = 1:length(fileNames)
    if fileNames(i).isdir == 0 % Check if it's not a directory
        fileName = fileNames(i).name;
        % Check if filename contains all strings in stringArray
        if all(contains(fileName, stringArray))
            break;
        end
    end
end
data_file_path = join([file_path,fileName],"");
data = readmatrix(data_file_path,"FileType","text","NumHeaderLines",1);
% distance 
r = data(:,1);
% Extract HTE0 fit from file (mu0's, T, lattice confinement and tunneling)
parameter_file_path = join([file_path,"Parameters_for_datasets.txt"],"");
[mu0, T, omega, tunneling] = find_mu0_vals(parameter_file_path,u,m);

%% Fitting data
xData = r;
yData = data(:,2:4);
initialGuess = [1,1,1,1];
order_max = 2;
otherParams = {t,u,m,n,dno,order_max,tunneling, omega};
[fittingParams, resnorm, residual, exitflag, output] = lsqcurvefit(@(fitParams,x) density_vs_r_fitting_function(x,fitParams,otherParams), initialGuess, xData, yData);
rmse = sqrt(mean(residual.^2));
%% mu-T grid
muq = zeros(numel(r),m);
for i=1:m
    muq(:,i) =  fittingParams(i+1) - .5*(6*1.66054e-27)*(2*pi*omega)^2*(r*752*10^-9).^2/(6.62607015e-34*tunneling); 
    % 752 nm is the lattice spacing distance.
end
%muq = linspace(-25,5,100); % balanced mus. can be changed for unbalanced case.
Tarray = [fittingParams(1)];
mu_file = join(['../data/csv_files/N=',num2str(m),'/muq.csv']);
T_file = join(['../data/csv_files/N=',num2str(m),'/T.csv']);
writematrix(muq,mu_file);
writematrix(Tarray,T_file);
clear("mu_file","T_file");


%% Generate the NLCE sum
order_max = 5;
orderlist = 1:order_max;
[density] = NLCE_sum(t,u,m,n,dno,{"density"},orderlist,muq,Tarray,false); % 4D double [m,mu,T,order]
%% Plot
% Necessary parameters
temperature_cut = T; % units of t.
orderlist = [1 4 5]; % NLCE order list to be plotted
obslist = {'density'}; % Obseravable list to be plotted.
Tarray = readmatrix(join(['../data/csv_files/N=',num2str(m),'/T.csv']));    %write a function to find the closes T value in the grid
[T,Tidx] = findClosestValue(Tarray,temperature_cut); % Find the T closest to the input and the corresponding slice index.
figure;
plot_order = 5;
for i=1:m
    data_plot(r,data(:,i+1),sprintf("n_%d",i),data(:,i+m+1));
    data_plot(r,density(i,:,1,plot_order),sprintf("NLCE %d n_%d",[plot_order,i]))
end


% try
%     for obs_idx = 1:length(obslist)
%         obs = obslist{obs_idx};
%         save=true;
%         imbalance_sun_plots(T,Tidx,t,u,m,n,dno,obs,orderlist,save);
%     end
% catch
%     % Define parameters not in the current workspace.
%     t = 1; u = [1 10 10]; m = 3; n = -1; dno = -1;
%     for obs_idx = 1:length(obslist)
%         obs = obslist{obs_idx};
%         save=true;
%         imbalance_sun_plots(T,Tidx,t,u,m,n,dno,obs,orderlist,save);
%     end
% end


%clear('all');


%% Functions 
function varargout = NLCE_sum(t,u,m,n,dno,observables,orderlist,muq,Tarray,save_files)
    outputs = {};
    itemlist = cellfun(@(x) x, observables, 'UniformOutput',false);
    itemlist{end + 1} = 'fileio';
    for i=1:numel(itemlist)
        if strcmp(itemlist{i},"density")
            NLCE_density = zeros(m,size(muq,1),numel(Tarray),numel(orderlist));
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
            side1 = [side1' 1:l]; side2 = [side2' 1:l]; %on-site correlator
            filename = join(['../data/mat_files/ED_mat_files/N=',num2str(m),'/' key],'');
            filename = join([filename," t=", double(t), "u=", double(u), "n=", n, "m=", m, "D=", dno, "ED.mat"],' '); % name of file that saves the output
            if any(strcmp(itemlist, "correlator"))
                itemlist{end + 1} = side1;
                itemlist{end + 1} = side2;
            end
            try
                disp(['Try to find ED data for graph ', num2str(key), '...']);
                % clear('output');
                load(filename); % try to load result of the ED to one given graph
                disp(timer_count);
                disp(['ED data for graph ', num2str(key), ' found.']);
            catch
                disp(['ED data for graph ', num2str(key), ' not found. Start to calculate ED.']);
    
                %[timer_count, spectra, nsigmapermute, testn, nimatrix, ni2matrix, domatrix, ninj] = ED_solver(g{k},t,u,n,m,dno,itemlist{:});
                [timer_count, spectra, nsigmapermute, testn, nimatrix] = ED_solver(g{k},t,u,n,m,dno,itemlist{:});
                disp(['ED calculation for graph ', num2str(key), ' finished.']);
                disp(timer_count);
            end
    
            disp(join([key, 'l=', l, "t=", double(t), "u=", double(u), "m=", m, "D=", dno, 'ED measurement start']));
            disp(join(["Measured item :", itemlist{:}]));
            %nimatrix = 1 + zeros(length(nimatrix))
    
            try
                nlce_path = join(["../data/mat_files/NLCE_mat_files/N=",num2str(m),"/"],'');
                nlce_filename = join([key," t=", double(t), "u=", double(u), "n=", n, "m=", m, "D=", dno, "LDA"],' ');
                nlce_filename = strrep(nlce_filename,'.','p');
                nlce = join([nlce_path,nlce_filename,".mat"],"");
                load(nlce);
                
                %WARNING : This deletes previous data. Need a better check.
                if ~isequal(mulist, muq, 'RelTol', 1e-6)
                    disp("mu-list do not match, deleting nlce graph data...");
                    delete(nlce);
                    load(nlce);
                end
                if ~isequal(Tlist, Tarray, 'RelTol', 1e-6)
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
                        densityaccum = zeros(1,m);
                        density_sq_accum = zeros(1,m+1);
                        for spin_idx = 1:m
                            rho_matrix = cellfun(@(x) x{spin_idx},nimatrix,'UniformOutput',false);
                            %rho2_matrix = cellfun(@(x) x{spin_idx},ni2matrix,'UniformOutput',false);
                            densityaccum(spin_idx) = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, rho_matrix);
                            %density_sq_accum(spin_idx) = LDA_measure(spectra,deltamu,T,nsigmapermute,testn, rho2_matrix);
                        end
                        % density_sq_accum(m+1) = LDA_measure(spectra,deltamu,T,nsigmapermute,testn,cellfun(@(x) x{m+1},ni2matrix,'UniformOutput',false));
                        % doaccum = zeros(1,uint8(m*(m-1)/2));
                        % % onsite_cor_accum = zeros(1,uint8(m*(m-1)/2));
                        % ninjaccum = zeros(1,uint8(m*(m-1)/2));
                        % for pair_idx = 1:uint8(m*(m-1)/2)
                        %     pair_matrix = cellfun(@(x) x{pair_idx},domatrix,'UniformOutput',false);
                        %     doaccum(pair_idx) = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, pair_matrix);
                        %     ninj_matrix = cellfun(@(x) x{pair_idx}, ninj,'UniformOutput',false);
                        %     correlatoraccum = (LDA_measure(spectra, deltamu, T, nsigmapermute, testn, ninj_matrix));
                        %     niniaccum(pair_idx) = sum(correlatoraccum(numEdge+1:end)); %on-site
                        %     if l>1
                        %         ninjaccum(pair_idx) = sum(correlatoraccum(1:numEdge)) - 1/l^2*numEdge*densityaccum(rev_pair_idx(pair_idx,1))*densityaccum(rev_pair_idx(pair_idx,2)); %nearest-neighbor
                        %     end
                        % end
                        %energyaccum = LDA_measure( spectra, deltamu, T, nsigmapermute, testn, 'energy') + mu0*sum(densityaccum) ;
                        %energy_sq = LDA_measure(spectra ,deltamu, T, nsigmapermute, testn, 'energy_square');
                        %cv_accum = 1/T^2*(energy_sq - energyaccum^2);
                        %entropyaccum = LDA_measure( spectra, deltamu, T, nsigmapermute, testn, 'entropy') ;
                        %cnnaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, cnn) ;
                        %ninjaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, ninj);
                        %  CTlist(end + 1, :) = [T mu0 densityaccum doaccum energyaccum entropyaccum cnnaccum' ninjaccum'];
                        alist(end + 1, :) = [densityaccum(:)'];
                    end
                end
                % nlce = join(["../data/mat_files/NLCE_mat_files/N=",num2str(m),"/" ,key],'');
                % nlce = join([nlce, " t=", double(t), "u=", double(u), "n=", n, "m=", m, "D=", dno, 'LDA.mat'],' ');
                mulist = muq; %To ensure the correct mus are used
                Tlist = Tarray;
                save(nlce,'alist','mulist','Tlist');
            end
            CTlist = CTlist + coefficient(k)*alist; %nlce sum
            disp(join([key, 'l=', l, "t=", double(t), "u=", double(u), "m=", m, "D=", dno, 'ED measurement finish']));
        end

        % Return the desired data.
        if exist("NLCE_density",'var') == 1
            for idx = 1:m
                NLCE_density(idx,:,:,order) = reshape(CTlist(:, idx),[size(muq,1),numel(Tarray)])';
            end
        end
        
        % Save data in files
        if save_files
            for idx = 1:m
                writematrix(reshape(CTlist(:, idx),[size(muq,1),numel(Tarray)])', join(['../data/csv_files/density_m'  num2str(m) '_flav' num2str(idx) '_u' num2str(double(u)) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']));
                % writematrix(1/T*reshape(CTlist(:, m+idx),[numel(muq),numel(Tarray)])', join(['../data/csv_files/compressibility_m'  num2str(m) '_flav' num2str(idx) '_u' num2str(double(u)) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']));
            end
            % writematrix(reshape(CTlist(:, 2*m + 1) - (sum(CTlist(:,1:m),2).^2),[numel(muq),numel(Tarray)])', join(['../data/csv_files/total_density_fluc_m'  num2str(m) '_u' num2str(double(u)) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']));
            % for idx = 1:uint8(m*(m-1)/2)
            %     connected = CTlist(:,rev_pair_idx(idx,1)).*CTlist(:,rev_pair_idx(idx,2));
            %     writematrix(reshape(CTlist(:,2*m+1+idx),[numel(muq),numel(Tarray)])', join(['../data/csv_files/doublon_m'  num2str(m) '_pair' num2str(idx) '_u' num2str(double(u)) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']));
            %     writematrix(reshape(CTlist(:,2*m+1+uint8(m*(m-1)/2)+3+idx)-connected,[numel(muq),numel(Tarray)])', join(['../data/csv_files/nini_m'  num2str(m) '_pair' num2str(idx) '_u' num2str(double(u)) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']));
            %     if order == 1
            %         connected = 0;
            %     end
            %     writematrix(reshape(CTlist(:,2*m+1+2*uint8(m*(m-1)/2)+3+idx),[numel(muq),numel(Tarray)])', join(['../data/csv_files/ninj_m'  num2str(m) '_pair' num2str(idx) '_u' num2str(double(u)) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']));
            % end
            % writematrix(reshape(CTlist(:,2*m+1+ uint8(m*(m-1)/2) + 2),[numel(muq),numel(Tarray)])', join(['../data/csv_files/entropy_m'  num2str(m) '_u' num2str(double(u)) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']));
            % writematrix(reshape(CTlist(:,2*m+1+ uint8(m*(m-1)/2) + 3),[numel(muq),numel(Tarray)])', join(['../data/csv_files/cv_m'  num2str(m) '_u' num2str(double(u)) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']));
        end
    end
    varargout = {NLCE_density};
end

function density = density_vs_r_fitting_function(x, fitParams, otherParams)
    % Assumes the ED files exist. Will not re-run the ED.
    T = fitParams(1);
    mu0 = fitParams(2:end);

    t = otherParams{1}; u = otherParams{2}; m = otherParams{3}; n = otherParams{4};
    dno = otherParams{5}; order_max = otherParams{6}; tunneling = otherParams{7}; omega = otherParams{8};
    
    mu = mu0 - .5*(6*1.66054e-27)*(2*pi*omega)^2*(x*752*10^-9).^2/(6.62607015e-34*tunneling); 
    disp(mu);

    density = zeros();
    for order=1:order_max
        [g, coefficient] = NLCE_load(order);
        for k = 1:numel(g)
            key = key_gen(g{k});
            l = numnodes(g{k});
            filename = join(['../data/mat_files/ED_mat_files/N=',num2str(m),'/', key],'');
            filename = join([filename," t=", double(t), "u=", double(u), "n=", n, "m=", m, "D=", dno, "ED.mat"],' '); % name of file that saves the output

            disp(['Try to find ED data for graph ', num2str(key), '...']);
            load(filename,'nimatrix','nsigmapermute','spectra','testn'); % try to load result of the ED to one given graph
            disp(['ED data for graph ',num2str(key), ' found.']);
            temp_list = [];
            for i = 1:numel(x)
                deltamu =  - mu(i,:); % REWRITE To include different mus
                densityaccum = zeros(1,m);
                for spin_idx = 1:m
                    rho_matrix = cellfun(@(x) x{spin_idx},nimatrix,'UniformOutput',false);
                    densityaccum(spin_idx) = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, rho_matrix);
                end
                temp_list(end + 1,:) = [densityaccum(:)'];
            end
            density = density + coefficient(k)*temp_list; %nlce sum
            disp(join([key, 'l=', l, "t=", double(t), "u=", double(u), "m=", m, "D=", dno, 'ED measurement finish']));
        end
    end
end

function imbalance_sun_plots(T,Tidx,t,u,m,n,dno,obs,orderlist,save)
    muq = readmatrix(join(['../data/csv_files/N=',num2str(m),'/muq.csv']));
    colorlist = {'g','c','y','r','b'};
    markerlist = {'o', 'x', '^'};
    rev_pair_idx = [1 2; 1 3; 2 3];

    % Setting up the y-range dictionary
    ylim_list = containers.Map;
    ylim_list('density') = [0,m];
    ylim_list('total_density_fluc') = [0,2];
    ylim_list('compressibility') = [0,1];
    ylim_list('doublon') = [0,m];
    ylim_list('entropy') = [0,2*log(m)];
    ylim_list('cv') = [-1,20];
    ylim_list('ninj') = [-.1,.1];
    ylim_list('nini') = [-.5,.5];

    % Setting up the y label dictionary.
    ylabel_list = containers.Map;
    ylabel_list('density') = '$\rho$';
    ylabel_list('doublon') = '$\mathcal{D}$';
    ylabel_list('total_density_fluc') = '$\langle\Delta \rho^2\rangle$';
    ylabel_list('compressibility') = '$\kappa$';
    ylabel_list('entropy') = '$\mathcal{S}$';
    ylabel_list('cv') = '$C_V$';
    ylabel_list('ninj') = '$\langle n_{a,i}n_{b,i+1}\rangle - \langle n_{a,i}\rangle\langle n_{b,i+1}\rangle$';
    ylabel_list('nini') = '$\langle n_{a,i}n_{b,i}\rangle - \langle n_{a,i}\rangle\langle n_{b,i}\rangle$';

    %Setting up the legend dictionary
    legend_list = containers.Map;
    legend_list('density') = 'flavor';
    legend_list('total_density_fluc') = '';
    legend_list('compressibility') = 'flavor';
    legend_list('doublon') = 'pair';
    legend_list('entropy') = '';
    legend_list('ninj') = 'pair';
    legend_list('nini') = 'pair';
    
    figure;
    for order=orderlist
        if strcmp(obs,'density') || strcmp(obs,'compressibility')
            len = m;
        elseif strcmp(obs,'doublon') || strcmp(obs,'ninj') ||strcmp(obs,'nini')
            len = uint8(m*(m-1)/2);
        else
            len = 1;
        end
        observable = zeros(length(muq),len);
        for idx = 1:len
            if strcmp(obs,'density') || strcmp(obs,'compressibility')
                fname = join(['../data/csv_files/' obs '_m'  num2str(m) '_flav' num2str(idx) '_u' num2str(double(u)) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']);
            elseif strcmp(obs,'doublon') || strcmp(obs,'ninj') || strcmp(obs,'nini')
                fname = join(['../data/csv_files/' obs '_m'  num2str(m) '_pair' num2str(idx) '_u' num2str(double(u)) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']);
            else
                fname = join(['../data/csv_files/' obs '_m'  num2str(m) '_u' num2str(double(u)) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']);
            end
            dataFromFile = readmatrix(fname);
            observable(:,idx) = dataFromFile(Tidx,:);
            if order==max(orderlist) && (strcmp(obs,'density') || strcmp(obs,'compressibility'))
                plot(muq,observable(:,idx),colorlist{order},'Marker',markerlist{idx},'MarkerSize',10,'LineStyle','none','DisplayName',join([legend_list(obs), num2str(idx)]));
            elseif order==max(orderlist) && (strcmp(obs,'doublon') || strcmp(obs,'ninj') || strcmp(obs,'nini'))
                plot(muq,observable(:,idx),colorlist{order},'Marker',markerlist{idx},'MarkerSize',10,'LineStyle','none','DisplayName',join([legend_list(obs), num2str(rev_pair_idx(idx,:))]));
            else
                h = plot(muq,observable(:,idx),colorlist{order},'Marker',markerlist(idx),'MarkerSize',10,'LineStyle','--','LineWidth',2);
                set(h,'HandleVisibility','off');
            end
            hold on;
        end
        if strcmp(obs,'density') || strcmp(obs,'doublon')
            if order==max(orderlist)
                plot(muq,sum(observable,2),'Color',colorlist{order},'LineStyle','--','LineWidth',2,'DisplayName','sum');
            else
                h = plot(muq,sum(observable,2),'Color',colorlist{order},'LineStyle','--','LineWidth',2);
                set(h,'HandleVisibility','off');
            end
        end
    end
    ylabel(ylabel_list(obs),'Interpreter','latex');
    if strcmp(obs,'entropy')
        plot([min(muq),max(muq)],[log(m),log(m)],'k-','LineWidth',2,'DisplayName',join(['log(',num2str(m),')']));
    end
    if strcmp(obs,'doublon')
        plot([min(muq),max(muq)],[1/m,1/m],'k-','LineWidth',2,'DisplayName','1/N');
    end
    grid on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    ax.Title.FontSize = 20;
    ylim(ylim_list(obs));
    xlim([min(muq),max(muq)]);
    xlabel('\mu');
    lgd = legend;
    lgd.FontSize = 16;
    lgd.Location = 'northwest';
    annotation('textbox', [0.35, 0.8, 0.1, 0.1], 'String', 'green: Atomic Limit', 'EdgeColor', 'none', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'green');
    annotation('textbox', [0.35, 0.75, 0.1, 0.1], 'String', 'red: NLCE 4', 'EdgeColor', 'none', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'red');
    annotation('textbox', [0.35, 0.7, 0.1, 0.1], 'String', 'blue: NLCE 5', 'EdgeColor', 'none', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'blue');
    title(join(["T/t=", num2str(T), ", U/t=(", num2str(u),")"]));
    ustr = num2str(u);
    ustr = strrep(ustr, '.', 'p');
    if save
        print(join(['../plots/',obs,'_u=',ustr,'_T=',num2str(T)]),'-dpng','-r300');
    end
end

function data_plot(xvals, yvals,varargin)
    if nargin<3
        yerr = [];
        legend_string = "";
    elseif nargin==3
        yerr = [];
        legend_string = varargin{1};
    else
        yerr = varargin{2};
        legend_string = varargin{1};
    end
    if isempty(yerr)
        plot(xvals,yvals,'x-','DisplayName',legend_string)
    else
        errorbar(xvals,yvals,yerr,'o-','DisplayName',legend_string)
    end
    grid on;
    legend('Location','northeast');
    hold on;
end

function [closestValue, index] = findClosestValue(array, target)
    % Calculate absolute differences between target and each element in array
    differences = abs(array - target);
    
    % Find the index of the element with the minimum difference
    [~, index] = min(differences);
    
    % Return the value from array corresponding to the index
    closestValue = array(index);
end

function [mu0, T, omega, tunneling] = find_mu0_vals(filepath, u, m)
    
    % Open the file for reading
    fileID = fopen(filepath, 'r');
    
    % Initialize variables to store the result
    mu0 = []; T = NaN; omega = NaN; tunneling = NaN;

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
    
    % Close the file
    fclose(fileID);
end

