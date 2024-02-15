[g,coefficients] = NLCE_load(5);
t=vpa(1); u=vpa(15.3); m=2; n=-1; dno=-1; mu=-u*(m-1)/2;
itemlist={'density','doubleoccupancy','fileio'};
[spectra, nsigmapermute, testn, nimatrix, domatrix] = ED_solver(g{1}, t, u, mu, n, m, dno,itemlist{:});
