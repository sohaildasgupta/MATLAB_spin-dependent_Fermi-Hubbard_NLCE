function varargout = ED_solver(graph, t, u, n, m, dno, varargin)
%
% ED_solver - solve ED of a given graph or lattice size with different parameters and arguments
%
% SYNTAX :n
%         [spectra, nsigmapermute] = ED_solver(graph, t, u, n, m, dno)
%         [spectra, nsigmapermute, eigenstates, basis] = ED_solver(graph, t, u, n, m, dno)
%         [spectra, nsigmapermute, testn, matrix] = ED_solver(graph, t, u, n, m, dno, observables)
%         [spectra, nsigmapermute, testn, matrix] = ED_solver(dim, t, u, n, m, dno, observables(if containing 'symmetries' option))
%         [gsenergy, groundstate, gsbasis] = ED_solver(dim, t, u, n, m, dno, observables(if containing 'groundstates' option)) - under construction
%
% INPUT ARGUMENTS :
% graph : input graph. If is a graph object, then the funciton reads he object; if is a lattice size vector [sz1, sz2, sz3], then the function automatically generates a graph of that dimension with periodci boundary condtion. NOTE: size vector dimension should equal to space dim, 1D 6x1 chain should bw wrtten as [6] instead of [6 1].
% t, U : parameters of FHM
% n : total particle #. n=-1 means no constraint on particle #.
% m : # of spin flavor species
% dno : truncation condition, when a state of the basis has states with onsite energy larger than dno*U, the basis state will be thrown away. dno=-1 means no truncation.
%
% observables : 'XXZ', 'density', 'doubleoccupancy', 'groundstates', 'fileio', 'symmetries', 'correlator'
% NOTE : 1. 'correlator' contains two inputs & two outputs
%            inputs : list of postion i & postion j
%            outputs : state expectation values of cnn & ninj
%        2. 'groundstates' calculates system ground state using Lanczos method, under construction
%        3. 'fileio' is the option to read results of diagonalization from and store them to the disk immediately instead of keeping them at the mem and store altogether after all sectors are calculated
%        4. 'symmetries' is the option to utilize spatial symmetries (so far only translation symmetries) when diagonalizing PBC rectangular system (or other lattice system in the future). When this option is on, the graph input switches to lattice dimension input, i.e.: graph object -> [lx ly lz] (array length eq to space dim)
%        5. 'XXZ' is the special mode to transfer the model calcuated to XXZ spin model according to the geometry equivalence between SU(2) hard-core Bose Hubbard model and XXZ spin model. Note that this transform is not quite efficient (only some simple param change, no special optimization). Under construction for output part.
%
% OUTPUT ARGUMENTS :
% spectra : energy spectrum of calculated system
% eigenstates : wavefunction of corresponding energy level
% basis : the basis used in the ED calculation
% nsigmapermute : information of spin permutation
% testn : total ptcl # of the sector, stored for calculating results changing mu & T without diagonalizing Hamiltonian again i.e.: thermodynamics calculation
% matrix : state expectation values of every eigenstates, stored for calculating results changing mu & T without diagonalizing Hamiltonian again i.e.: thermodynamics calculation
%
% NOTE : 1. The basis states are ordered by spin flavor first, then by site index. e.g.: if |ud> = C^dagger_1u C^dagger_2d|0>, then |du> = - C^dagger_2u C^dagger_1d|0>, so the SU(2) singlet |ud>-|du>=C^dagger_1u C^dagger_2d|0> + C^dagger_2u C^dagger_1d|0>. As a result, the SU(2) spin singlet has the same phase factor for both components.
% 2. For memory saving please use 'fileio' to directly use disk space for saving output step-by-step. Remember that the full eigenstate matrix is the most memory-consumming output. It has the same size as the Hamiltonian matrix while it is not sparse.
%
% Use graph key for NLCE summation.
% Author : HT Wei, 2019 @ Rice
%
[ni, ni2, do, p2, nn, io, sun_symm] = deal(false);
itemlist = {};
if all(u == u(1))
    sun_symm = true;
end
sun_symm = false;
timer_count = tic;

try
    l = numnodes(graph);
catch
    latticedims = graph;
    spacedim = numel(latticedims);
    period = latticedims;
    graph = gridgraph(graph, 'pbc');
    l = numnodes(graph);
end

if numel(varargin)
    input = varargin;
    for idx = 1:numel(input)
        arg = input{idx};
        if strcmp(arg, 'density')
            ni = true;
            itemlist{end + 1} = arg;
            continue
        elseif strcmp(arg, 'density2')
            ni2 = true;
            itemlist{end + 1} = arg;
            continue
        elseif strcmp(arg, 'doubleoccupancy')
            do = true;
            itemlist{end + 1} = arg;
            continue    
        elseif strcmp(arg, 'fileio')
            io = true;
            itemlist{end + 1} = arg;
            continue
        elseif strcmp(arg, 'correlator')
            nn = true;
            itemlist{end + 1} = arg;
            posi = input{idx + 1};
            posj = input{idx + 2};
            % poslen = length(posi);
            break
        end
    end
end

% sectors = int8(-1); % To cancel the for loop of different sectors
% symopt = int8(0);

key = key_gen(graph);
disp(join([key, "l=", l, "t=", double(t), "u=", double(u), "n=", n, "m=", m, "D=", dno, "ED start"]));
[side1, side2] = findedge(graph);
orderlist = [side1 side2];

disp(join(["Output options :", itemlist]));
filename = join(['../data/mat_files/ED_mat_files/N=',num2str(m),'/' key],'');
filename = join([filename," t=", double(t), "u=", double(u), "n=", n, "m=", m, "D=", dno, 'ED.mat'],' '); % name of file that saves the output

nsigmatuples = unique(sort(permn(0:l, m), 2), 'row'); % for
if dno >= 0
    if 0 <= n && n <= l + dno
        % nsigmatuples = nsigmatuples(sum(nsigmatuples, 2) == n, :);
        nsigmatuples = nsigmatuples(sum(nsigmatuples, 2) <= n, :);
    else
        nsigmatuples = nsigmatuples(sum(nsigmatuples, 2) <= l + dno, :);
    end
elseif 0 <= n && n <= l * m
    % nsigmatuples = nsigmatuples(sum(nsigmatuples, 2) == n, :);
    nsigmatuples = nsigmatuples(sum(nsigmatuples, 2) <= n, :);
end

% Generate all the spin permutations
nsigma = {}; 
total_spin_permute_no = 0; %Counts the total number of spin permutations
for spinnocofig = nsigmatuples' %matlab iterates over columns
    spin_permutations = unique(perms(spinnocofig), 'row'); 
    total_spin_permute_no = total_spin_permute_no + size(spin_permutations,1);
    nsigma{end + 1} = spin_permutations;
end
if sun_symm
    total_spin_permute_no = length(nsigma);
end

nsigmalen = length(nsigma);
sitenotuples = permn(logical(0:1), l);
subbasisvectorlist = cell(1, l + 1);

for particleno = 0:l
    subbasisvectorlist{particleno + 1} = sitenotuples(sum(sitenotuples, 2) == particleno, :);
end

clear('sitenotuples');

    function hfunction()
        %
        % Generate tuneling matrix for given flavor and particle number
        %
        subbasisvector = cell2mat(subbasisvector);
        dim = size(subbasisvector, 1);
        idx = []; jdx = []; value = [];

        for cidx = 1:dim
            startstate = subbasisvector(cidx, :);
            for swap = orderlist'
                if sum(startstate(swap)) ~= 1
                    continue
                end
                endstate = startstate;
                endstate(swap') = startstate(flip(swap'));
                ridx = find(all(subbasisvector(1:cidx - 1, :) == endstate, 2));
                if ridx < cidx% Only construct upper triangular part
                    connection = sort(swap);
                    % To give the right sign for correction indexing of basis states
                    v = -(-1)^sum(startstate((connection(1) + 1):(connection(end) - 1))) * t;
                    [idx(end + 1), jdx(end + 1), value(end + 1)] = deal(ridx, cidx, v);
                    % [idx(end + 1), jdx(end + 1), value(end + 1)] = deal(cidx, ridx, v);
                end
            end
        end
        hoppinglist{end + 1} = sparse(idx, jdx, value, dim, dim);
    end

hoppinglist = {};
for subbasisvector = subbasisvectorlist
    hfunction();
end

spectra = {}; nsigmapermute = []; % time = 0;
nOutput = nargout;

if nOutput > 2
    if numel(varargin)
        testn = {};
        if ni
            nimatrix = {};
        end
        if ni2
            ni2matrix = {};
        end
        if do
            domatrix = {};
        end
        if nn
            cnn = {};
            ninj = {};
        end
        if io
            save(filename, 'spectra', 'nsigmapermute', '-v7.3');
        end
    else
        eigenstates = {}; basis = {};
    end
end
    function varmem()
        %
        % Show memory usage of the 5 largest variables.
        %
        bytecount = whos();
        varname = {bytecount.name};
        [bytecount, byteorder] = maxk([bytecount.bytes] / 1048576, 5);
        disp(join(["Memory used by variables [", varname{byteorder}, "] = [", bytecount, "] MBytes"]));
    end


permute_counter=0;
for scidx = 1:nsigmalen
    % tic
    spin_configs = nsigma{nsigmalen + 1 - scidx};
    if sun_symm
        permuteno = size(spin_configs,1);
        spinconfig_len = 1;
    else
        permuteno = 1;
        spinconfig_len = size(spin_configs,1);
    end
    
    for spin_config_idx = 1:spinconfig_len
        permute_counter = permute_counter+1;
        spinnocofig = spin_configs(spin_config_idx, :);
        disp(join([permute_counter, "/", total_spin_permute_no, "nsigma species = [", spinnocofig, "] start diagonalizing"]));
        basisvector = cell(1, m); hmatrix = cell(1, m); dimlist = ones(1, m);

        for spidx = 1:m
            basisvector{spidx} = subbasisvectorlist{spinnocofig(spidx) + 1};
            hmatrix{spidx} = hoppinglist{spinnocofig(spidx) + 1};
            dimlist(spidx) = size(basisvector{spidx}, 1);
        end

        fockdim = prod(dimlist);
        basisvector = reshape(celltuples(basisvector), fockdim, l, m);
        % Kron module for tunneling matrix
        if m > 1
            tmatrix = sparse(fockdim, fockdim);
            for idx = 1:m
                if idx == 1
                    C = hmatrix{1};
                    for jdx = 2:m
                        C = kron(C, speye(dimlist(jdx)));
                    end
                else
                    C = speye(dimlist(1));
                    for jdx = 2:m
                        if jdx == idx
                            C = kron(C, hmatrix{jdx});
                        else
                            C = kron(C, speye(dimlist(jdx)));
                        end
                    end
                end
                tmatrix = tmatrix + C;
            end
        else
            tmatrix = hmatrix{1};
        end

        clear('hmatrix', 'C');

        onsiteparticleno = sum(basisvector,3);
        %module for the interaction matrix and particle numbers
        if ~sun_symm
            onsiteinteraction = zeros(fockdim,1);
            pair_count = zeros(uint8(m*(m-1)/2),fockdim);
            for fockdim_idx = 1:fockdim
                for site_idx = 1:l
                    flavor_index = find(basisvector(fockdim_idx,site_idx,:));
                    if numel(flavor_index)>1
                        index_pairs = nchoosek(flavor_index,2);
                        for pair = index_pairs'
                            small = min(pair);
                            big = max(pair);
                            uidx = (small-1)*m - uint8(small*(small-1)/2) + big-small;
                            onsiteinteraction(fockdim_idx) = onsiteinteraction(fockdim_idx) + u(uidx);
                            pair_count(uidx,fockdim_idx) = pair_count(uidx,fockdim_idx) + 1;
                        end
                    end
                end
            end
        else
            onsiteinteraction = u(1) * dot(onsiteparticleno, onsiteparticleno - 1, 2) / 2;
        end
        p2_op = eq(onsiteparticleno,2);

        if dno >= 0
            truncatepos = find(onsiteinteraction <= dno);
            onsiteparticleno = onsiteparticleno(truncatepos, :);
            onsiteinteraction = onsiteinteraction(truncatepos);
            basisvector = basisvector(truncatepos, :, :);
            truncatedfockdim = length(truncatepos);
            if l > 1
                tmatrix = tmatrix(truncatepos, truncatepos);
            end
        else
            truncatedfockdim = fockdim;
        end
        % for k = sectors'
        %     k = double(k);
            ham = tmatrix;
            clear('tmatrix');
            if ~truncatedfockdim
                disp('Hamiltionian matrix empty, skip to the next sector.');
                continue
            end

            % NOTE : the matrix is only upper triangular
            ham = ham + ham'; % ' is already Hermitian conjugate
            ham = ham + spdiags(onsiteinteraction, 0, truncatedfockdim, truncatedfockdim);% - (mu + u * (m-1) / 2) * sum(spinnocofig) * speye(truncatedfockdim);
            ham = full(ham);
         
            disp(join(["Matrix dim =", truncatedfockdim]));
            varmem();

            disp('Start diagonalizing Hamiltonian matrix...');
            % tic
            if nOutput > 2
                [nsigmaeigenstates, nsigmaspectra] = eig(ham, 'vector'); % NOTE : nsigmaeigenstates are column vectors
                % [nsigmaeigenstates, nsigmaspectra] = eigs(ham, truncatedfockdim);
            else
                nsigmaspectra = eig(ham);
            end
            % toc
            clear('ham');

            % Load file in order to append new results to output variables
            if io
                try
                    load(filename);
                end
            end
            if numel(varargin)
                testn{end + 1} = spinnocofig;
                if nn
                     % if l == 1
                    %     cnn{end + 1} = zeros(1, truncatedfockdim);
                    %     % ninj{end + 1} = zeros(1, truncatedfockdim);
                     % else
                        
                        if ~sun_symm
                            ninj_mean = {};
                            for idx = 1:m
                                for jdx = idx+1:m
                                    ninj_mean{end + 1} = (basisvector(:,posi,idx) .* basisvector(:,posj,jdx))' * abs(nsigmaeigenstates).^2;
                                end
                            end
                            ninj{end + 1} = ninj_mean;
                            clear('ninj_mean');
                        else
                            cnnmatrix = 0;
                            for spinidx = 1:m
                                spina = repmat(spinidx, [1 m-1]);
                                spinb = [1:spinidx-1 spinidx+1:m];
                                cnnmatrix = cnnmatrix + sum(basisvector(:, posi, spina) .* (basisvector(:, posj, spina) - basisvector(:, posj, spinb)), 3);
                            end
                            ninjmatrix = permuteno * (onsiteparticleno(:, posi) .* onsiteparticleno(:, posj))' * abs(nsigmaeigenstates).^2;
                            cnnmatrix = permuteno * cnnmatrix' * abs(nsigmaeigenstates).^2;
                            cnn{end + 1} = cnnmatrix;
                            ninj{end + 1} = ninjmatrix;
                            clear('cnnmatrix', 'ninjmatrix');
                        end
                        
                    
                end
                if ni
                    if ~sun_symm
                        ni_mean = {};
                        for spin_idx = 1:m
                            ni_mean{end+1} = sum(basisvector(:,:,spin_idx),2)' * abs(nsigmaeigenstates).^2;
                        end
                        nimatrix{end + 1} = ni_mean;
                    else
                        nimatrix{end+1} = permuteno * sum(onsiteparticleno,2)' * abs(nsigmaeigenstates).^2;
                    end
                end
                if ni2 %assuming no sun_symm
                    ni2_mean = {};
                    for spin_idx = 1:m
                        ni2_mean{end+1} = sum(basisvector(:,:,spin_idx).^2,2)' * abs(nsigmaeigenstates).^2;
                    end
                    ni2_mean{end + 1} = sum((onsiteparticleno).^2,2)' * abs(nsigmaeigenstates).^2;
                    ni2matrix{end + 1} = ni2_mean;
                end
                if do
                    if ~sun_symm
                        do_mean = {};
                        for pair_idx = 1: uint8(m*(m-1)/2)
                            do_mean{end+1} = pair_count(pair_idx,:) * abs(nsigmaeigenstates).^2;
                        end
                        domatrix{end + 1} = do_mean;
                    else
                        domatrix{end + 1} = permuteno * onsiteinteraction' * abs(nsigmaeigenstates).^2;
                    end
                end
            elseif nOutput > 2
                eigenstates{end + 1} = nsigmaeigenstates;
                basis{end + 1} = basisvector;
            end

            spectra{end + 1} = nsigmaspectra;
            nsigmapermute(end + 1) = permuteno;

            clear('nsigmaspectra', 'nsigmaeigenstates');
            if io
                save(filename, 'spectra', 'nsigmapermute', 'testn', '-v7.3');

                if ni
                    save(filename, 'nimatrix', '-append');
                end
                if ni2
                    save(filename, 'ni2matrix', '-append');
                end
                if do
                    save(filename, 'domatrix', '-append');
                end

                if nn
                    save(filename, 'cnn', 'ninj', '-append');
                end

                % if scidx ~= nsigmalen
                clear('spectra','nsigmapermute', 'testn', 'nimatrix', 'domatrix', 'cnn', 'ninj');
                % end
            end
        % end
        clear('tmatrix', 'basisvector', 'onsiteparticleno', 'originalonsitptclno', 'onsiteinteraction');
    end
end

if io
    try
        load(filename);
    end
end
timer_count = toc(timer_count);
outputs = {timer_count};
outputs{end+1} = spectra;
outputs{end + 1} = nsigmapermute;

if nOutput > 2
    if numel(varargin)
        outputs{end + 1} = testn;
        if ni
            outputs{end + 1} = nimatrix;
        end
        if ni2
            outputs{end + 1} = ni2matrix;
        end
        if do
            outputs{end + 1} = domatrix;
        end
        if nn
            % outputs{end + 1} = cnn;
            outputs{end + 1} = ninj;
        end
    else
        outputs = [outputs, {eigenstates}, {basis}];
    end
end
if ~io
    save(filename, 'outputs', '-v7.3'); % save all result of the ED
end

varargout = outputs;

disp(join([key, "l=", l, "t=", double(t), "u=", double(u), "n=", n, "m=", m, "D=", dno, "ED finish"]));
end

function counter = sortcount(C)
%
% Count # of n.n. Fermion swaps when ordering state.
%
counter = 0;

for index = 1:size(C, 3)
    exchangelist = nonzeros(C(:, :, index))';
    if length(exchangelist) > 1
        swaped = true;
        while swaped
            swaped = false;
            for idx = 1:length(exchangelist) - 1
                if exchangelist(idx) > exchangelist(idx + 1)
                    exchangelist([idx idx+1]) = exchangelist([idx+1 idx]);
                    counter = counter + 1;
                    swaped = true;
                end
            end
        end
    end
end
end

function sgn = symmsgn(symbsvector, symopt)
%
% Determine the sign of the state vector ater symmetry operation
%
sgn = symbsvector .* (1:size(symbsvector, 2));
sgn = sgn(:, symopt, :);
sgn = (-1)^sortcount(sgn);
end

function [trubsvec, onsiteparticleno, onsiteinteraction] = par_vec_select(testbsvec, lastbsvec, dimlist, m, dno, symm, symopt)
nsym = size(symopt, 1);
trubsvec = logical.empty; onsiteparticleno = []; onsiteinteraction = [];
idxlist = ones(1, m - 1);
ready = false;
if ~isempty(lastbsvec)
    while ~ready
        for spidx = 2:m
            testbsvec(1, :, spidx) = lastbsvec{spidx - 1}(idxlist(spidx - 1), :);
        end
        
        % Update index:
        ready = true; % Assume that the WHILE loop is ready
        for k = 1:m - 1
            idxlist(k) = idxlist(k) + 1;
            if idxlist(k) <= dimlist(k + 1)
                ready = false;
                break% v(k) increased successfully, leave "for k" loop
            end
            idxlist(k) = 1; % v(k) reached the limit, reset it
        end

        % Select states
        testptclno = sum(testbsvec, 3);
        testinteraction = testptclno * (testptclno - 1)' / 2;
        
        if dno >= 0 && testinteraction > dno
            continue
        end
        
        trubsvec(end + 1, :, :) = testbsvec;
        onsiteparticleno(end + 1, :) = testptclno;
        onsiteinteraction(end + 1, 1) = testinteraction;
    end
end
end
