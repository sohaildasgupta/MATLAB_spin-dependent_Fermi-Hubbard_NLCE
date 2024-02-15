feature('numcores');

gridsize = 6;
l = prod(gridsize);
graph = gridgraph(gridsize, 'obc');
key = key_gen(graph);
t = 1; u = 2;
m = 2;
n = -1;
dno = -1;

dT = 0.1;
Tmax = 3;
n0 = 2.4e+4; % Total particle # in experiment.
% Tq0 = readmatrix(['T_vals.csv'])';
% muq = readmatrix(['mu_vals.csv'])';
Tq0 = logspace(-1,2,100);


j0 = [2,3,4,5,6];
% for m = [2 6]
[side1, side2] = deal([1,1,1,1,1],[2,3,4,5,6]); % every n.n. correlation for PBC 1D chain
itemlist = {'density', 'doubleoccupancy', 'fileio', 'correlator', side1, side2};
CTlist = [];
muq = 1.;

mu = -u * (m - 1) / 2;


% for dno = 1
% for n = 4:6

filename = join([key, "t=", t, "u=", u, "mu=", mu, "n=", n, "m=", m, "D=", dno, 'ED.mat']);

try
    disp(['Try to find ED data for graph ', num2str(key), '...']);
    clear('output');
    load(filename); % try to load result of the ED to one given graph
    disp(['ED data for graph ', num2str(key), ' found.']);
catch
    disp(['ED data for graph ', num2str(key), ' not found. Start to calculate ED.']);
    tic
    [spectra, nsigmapermute, testn, nimatrix, domatrix, cnn, ninj] = ED_solver(graph, t, u, mu, n, m, dno, itemlist{:});
    toc
    disp(['ED calculation for graph ', num2str(key), ' finished.']);
end

tic
disp(join([key, 'l=', l, "t=", t, "u=", u, "m=", m, "D=", dno, 'ED measurement start']));
disp(join(["Measured item :", itemlist{:}]));

    
for T = Tq0

    for mu0 = muq
        deltamu = mu - mu0;
        densityaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, nimatrix); % particle # per site, needed to transform to per unit cell
        doaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, domatrix)/l;
        energyaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, 'energy') / l + (mu0 + u * (m - 1) / 2) * densityaccum; % entropy per site, needed to transform to per unit cell
        entropyaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, 'entropy') /l; % entropy per site, needed to transform to per unit cell
        cnnaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, cnn) ;
        ninjaccum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, ninj);
        CTlist(end + 1, :) = [T mu0 densityaccum doaccum energyaccum entropyaccum cnnaccum' ninjaccum'];
    end

end

disp(join([key, 'l=', l, "t=", t, "u=", u, "m=", m, "D=", dno, 'ED measurement finish']));
toc

subplot(2,1,1,'XScale','log'); hold on
grid on;
for j=j0
    txt = join(['C_{NN}(r=' num2str(j-1) ')'],'');
    semilogx(Tq0,CTlist(:,7+(j-2)),'.','DisplayName',txt); hold on;
end
legend show;

subplot(2,1,2,'XScale','log'); hold on
grid on;
for j=j0
    txt = join(['NiNj(r=' num2str(j-1) ')'],'');
    semilogx(Tq0,CTlist(:,12+(j-2)),'.','DisplayName',txt); hold on;
end
legend show;
%{
writematrix(CTlist, join([key, 'l=', l, "u=", u, "n=", n, "m=", m, "D=", dno, "Homo.xlsx"]));
writematrix(muq', 'mu_vals.csv');
writematrix(Tq0', 'T_vals.csv');
writematrix(reshape(CTlist(:, 3), length(muq), length(Tq0))', 'Densities.csv');
writematrix(reshape(CTlist(:, 4), length(muq), length(Tq0))', 'Doubleoccupancies.csv');
writematrix(reshape(CTlist(:, 5), length(muq), length(Tq0))', 'Energies.csv');
writematrix(reshape(CTlist(:, 6), length(muq), length(Tq0))', 'Entropies.csv');
%}
%writematrix(reshape(CTlist(:, 7), length(muq), length(Tq0))', 'CNNs.csv');
%writematrix(reshape(CTlist(:, end + 1 - length(side1)), length(muq), length(Tq0))', 'ninis.csv');
% adiabatic_loading
