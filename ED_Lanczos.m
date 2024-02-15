feature('numcores');

try
    pc = parcluster
    % parpool('local', str2num(getenv('SLURM_TASKS_PER_NODE')));
    parpool('local', feature('numcores'));
end

gridsize = [3 3];
l = prod(gridsize);
graph = gridgraph(gridsize, 'pbc');
% key = key_gen(graph);
key = gridsize;
graph = gridsize;
t = 1; u = 15.3;
m = 3;
n = l;
dno = 1;
mu = -u * (m - 2) / 2;
% M = 173 * 1.66054e-27;
% omega = 2 * pi^2 * 155 * 43 * 126; % 1D parameters needed to change, but now for comparison with 2D we use 3D parameters
% ht = 161 * 6.62607004e-34;
% L = 266e-9;
% constant = (ht / M)^(3/2) / (L^3 * omega);
% constant = 1169.94090666726; % 1D
% constant = 725.327812187453; % 3D

% dT = 0.02;
% dr = 0.05;
% maxr = 20; % T=15, r=15 there is about 4 added to ptcl #, negligible compared to its order of 2.4e4
% n0 = 2.4e+4; % Total particle # in experiment.
% Tq0 = dT:dT:3; % NOTE : the difference between Tq in this .m file and the ED_LDA_muT_table_gen.m file !!!
% [side1, side2] = deal([1 1]', [2 4]'); % every n.n. correlation for PBC 1D chain
itemlist = {'groundstates', 'fileio', 'symmetries'};

% for m = [2 6]
% CTlist = zeros(length(Tq), 4 + length(side1) * 2, 2);
CTlist = [];

% for dno = 1:2

%     for n = gridsize:gridsize + dno

filename = join(["l=", key, "t=", t, "u=", u, "mu=", mu, "n=", n, "m=", m, "D=", dno, 'GS.mat']);
% muTtablefilename = join(['1D', key, 'l=', l, 'm=', m, 'u=', u, "n=", n, "dno=", dno, 'mu-T table.xlsx']);

% try % find out data of interpolation
%     ipldata = readmatrix(muTtablefilename);
% catch
%     ED_LDA_muT_table_gen% NOTE : the difference between Tq in this .m file and the ED_LDA_muT_table_gen.m file !!!
%     writematrix(ipldata, muTtablefilename);
% end

try
    disp(['Try to find ED data for graph ', num2str(key), '...']);
    clear('output');
    load(filename); % try to load result of the ED to one given graph
    disp(['ED data for graph ', num2str(key), ' found.']);
catch
    disp(['ED data for graph ', num2str(key), ' not found. Start to calculate ED.']);
    tic
    [spectra, gs, gsb, k] = ED_solver(graph, t, u, mu, n, m, dno, itemlist{:});
    toc
    disp(['ED calculation for graph ', num2str(key), ' finished.']);
end

% STlist = [];
% f = fit(ipldata(:, 2), ipldata(:, 1), 'poly3');

% tic
% disp(join([key, 'l=', l, "t=", t, "u=", u, "m=", m, "D=", dno, 'ED LDA measurement start']));
% disp(join(["Measured item :", itemlist{:}]));

% for T = Tq0
%     mu0 = f(T);
%     ptclno = 0;
%     S = 0;
%     C_NN = 0;
%     NN = 0;

% for r = 0:dr:maxr
%         mu_prime = mu0 - r^2/2;
%         deltamu = mu - mu_prime;
%         densityaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, nimatrix); % particle # per site, needed to transform to per unit cell
%         entropyaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, 'entropy') / l; % entropy per site, needed to transform to per unit cell
%         cnnaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, cnn);
%         ninjaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, ninj);

%         % if r == 0
%         %     trapcenterdensity = densityaccum;
%         % end

%         if densityaccum > 1e-6
%             ptclno = ptclno + constant * densityaccum * r^2 * dr;
%             S = S + constant * entropyaccum * r^2 * dr;
%             C_NN = C_NN + constant * cnnaccum * r^2 * dr;
%             NN = NN + constant * ninjaccum * r^2 * dr;
%         else
%             break
%         end

%     end

%     A = -C_NN' / (2 * ptclno);
%     I = -2 * C_NN' ./ (NN' - C_NN');

%     if ~isnan(A)
%         STlist(end + 1, :) = [T mu0 ptclno S/ptclno A I];
%     end

% end

% disp(join([key, 'l=', l, "t=", t, "u=", u, "m=", m, "D=", dno, 'ED LDA measurement finish']));
% toc
% writematrix(STlist, join([key, 'l=', l, "u=", u, "n=", n, "m=", m, "D=", dno, "STO.xlsx"]));
% CTlist(:, :, end + 1) = STlist;
% % end

% writematrix(CTlist, join([key, 'l=', l, "u=", u, "STO.xlsx"]), 'Sheet', m / 2);
% % end
