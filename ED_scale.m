feature('numcores');
cores = [1 2 4 6];
Efficiency = [];
Speedup = [];
for i=1:numel(cores);
    maxNumCompThreads(cores(i));
    tic
    [g,coefficients] = NLCE_load(5);
    t=vpa(1); u=vpa(15.3); m=6; n=6; dno=3; mu=-u*(m-1)/2;
    itemlist={'density','doubleoccupancy','fileio'};
    [spectra, nsigmapermute, testn, nimatrix, domatrix] = ED_solver(g{1}, t, u, mu, n, m, dno,itemlist{:});
    walltime(i) = toc;
    Speedup(end+1) = walltime(1)/walltime(i);
    Efficiency(end+1) = 100*Speedup(i)/i;
end
%x = figure('visible','off');
%lologlog(cores,Speedup)
%exportgraphics(x,'scale.png')
newdata = [cores',Speedup'];
csvwrite('test.csv',newdata)