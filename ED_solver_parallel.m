function output = ED_solver_parallel(graph, t, u, n, m, dno, filename ,varargin)
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
% observables : 'density', 'doublon', 'triplon', 'onsite_pair', 'nearest_pair'
% NOTE : 1. 'correlator' contains two inputs & two outputs
%            inputs : list of postion i & postion j
%            outputs : state expectation values of cnn & ninj
%        2. 'groundstates' calculates system ground state using Lanczos method, under construction
%        3. 'fileio' is the option to read results of diagonalization from and store them to the disk immediately instead of keeping them at the mem and store altogether after all sectors are calculated
%        4. 'symmetries' is the option to utilize spatial symmetries (so far only translation symmetries) when diagonalizing PBC rectangular system (or other lattice system in the future). When this option is on, the graph input switches to lattice dimension input, i.e.: graph object -> [lx ly lz] (array length eq to space dim)
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
% Co-author : S Dasgupta, 2024 @ Rice
%
[ni, ni2, do, tri, p2, nn, sun_symm] = deal(false);
itemlist = {};
if isscalar(u)
    sun_symm = true;
end

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
        elseif strcmp(arg, 'doublon')
            do = true;
            itemlist{end + 1} = arg;
            continue   
        elseif strcmp(arg, 'triplon')
            tri = true;
            itemlist{end + 1} = arg;
            continue
        elseif strcmp(arg, 'fileio')
            io = true;
            itemlist{end + 1} = arg;
            continue
        elseif strcmp(arg, 'nearest_pair')
            nn = true;
            itemlist{end + 1} = arg;
            if l>1
                posi = input{idx + 1};
                posj = input{idx + 2};
            end
            % poslen = length(posi);
            continue
        elseif strcmp(arg, 'p2')
            p2 = true;
            itemlist{end + 1} = arg;
            continue
        end
    end
end

key = key_gen(graph);
disp(join([key, "l=", l, "t=", double(t), "u=", double(u), "n=", n, "m=", m, "D=", dno, "ED start"]));
[side1, side2] = findedge(graph);
orderlist = [side1 side2];

disp(join(["Output options :", itemlist]));

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
permuteno = [];

for spinnocofig = nsigmatuples' %matlab iterates over columns
    spin_permutations = unique(perms(spinnocofig), 'row'); 
    if sun_symm
        nsigma{end + 1} = spin_permutations(1,:); % Choose only the first permutation.
        permuteno(end + 1) = size(spin_permutations,1);
    else
        for permute_idx = 1:size(spin_permutations,1)
            nsigma{end + 1} = spin_permutations(permute_idx,:);
            permuteno(end + 1) = 1;
        end
    end
end
% Total number of independent spin-configurations.
total_spin_permute_no = length(nsigma); 

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
    i_dx = []; j_dx = []; value = [];

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
                [i_dx(end + 1), j_dx(end + 1), value(end + 1)] = deal(ridx, cidx, v);
                % [idx(end + 1), jdx(end + 1), value(end + 1)] = deal(cidx, ridx, v);
            end
        end
    end
    hoppinglist{end + 1} = sparse(i_dx, j_dx, value, dim, dim);
end

hoppinglist = {};
for subbasisvector = subbasisvectorlist
    hfunction();
end

spectra = cell(1,total_spin_permute_no); 

if numel(varargin)
    if sun_symm
        testn = zeros(1,total_spin_permute_no);
    else
        testn = cell(1,total_spin_permute_no);
    end
    if ni
        nimatrix = cell(1,total_spin_permute_no);
    end
    if ni2
        ni2matrix = cell(1,total_spin_permute_no);
    end
    if do
        domatrix = cell(1,total_spin_permute_no);
    end
    if tri
        trimatrix = cell(1,total_spin_permute_no);
    end
    if nn
        cnn = cell(1,total_spin_permute_no);
        ninj = cell(1,total_spin_permute_no);
    end
    if p2
        p2matrix = cell(1,total_spin_permute_no);
    end
    % if io
    % save(filename, 'spectra', 'nsigmapermute', '-v7.3');
    % end
else
    eigenstates = {}; basis = {};
end

% Compute the dimension of each sector
fockdimlist = ones(1,total_spin_permute_no);
dimlist_arr = ones(total_spin_permute_no,m);
for scidx = 1:total_spin_permute_no  
    basisvector = cell(1,m);
    spinnocofig = nsigma{scidx};
    for spidx = 1:m
        basisvector{spidx} = subbasisvectorlist{spinnocofig(spidx) + 1};
        dimlist_arr(scidx,spidx) = size(basisvector{spidx},1);
    end
    fockdimlist(scidx) = prod(dimlist_arr(scidx,:));
end

% Sort in ascending order of Fock-space dimension.
[fockdimlist_sort, fock_ind] = sort(fockdimlist);
nsigma_sort = nsigma(fock_ind);
dimlist_arr_sort = dimlist_arr(fock_ind,:);

%Find the index for the threshold
size_threshold = 5000; % Beyond this parallelize
threshold_index = find(fockdimlist_sort>size_threshold,1,'first');


% Run the small matrices serially.
%{
for scidx = 1:threshold_index-1
    spinnocofig = nsigma_sort{scidx};
    disp(join([scidx, "/", total_spin_permute_no, "nsigma species = [", spinnocofig, "] start diagonalizing"]));
    basisvector = cell(1, m); hmatrix = cell(1, m); 
    dimlist = dimlist_arr_sort(scidx,:);
    fockdim = fockdimlist_sort(scidx);
    for spidx = 1:m
        basisvector{spidx} = subbasisvectorlist{spinnocofig(spidx) + 1};
        hmatrix{spidx} = hoppinglist{spinnocofig(spidx) + 1};
    end 
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
        triplon_count = zeros(uint8(m*(m-1)*(m-2)/6),fockdim);
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
                if numel(flavor_index)>2 % works for only 3 flavors. update for general n-flavors.
                    triplon_count(1,fockdim_idx) = triplon_count(1,fockdim_idx) + 1;
                end
            end
        end
        truncatedfockdim = fockdim;
        % Include truncation
        %
        %
    else
        onsiteinteraction = dot(onsiteparticleno, onsiteparticleno - 1, 2) / 2;
               
        % Basis-state truncation
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
        % p2_op = eq(onsiteparticleno,2);
        triplon_count = sum(onsiteparticleno .*(onsiteparticleno-1) .* (onsiteparticleno-2),2)/6;
    end
  
    ham = tmatrix;
    clear('tmatrix');
    
    if ~truncatedfockdim
        disp('Hamiltionian matrix empty, skip to the next sector.');
        continue
    end

    % NOTE : the matrix is only upper triangular
    ham = ham + ham'; % ' is already Hermitian conjugate
    if sun_symm
        ham = ham + u*spdiags(onsiteinteraction, 0, truncatedfockdim, truncatedfockdim);
    else
        ham = ham + spdiags(onsiteinteraction, 0, truncatedfockdim, truncatedfockdim);
    end
    ham = full(ham);
 
    disp(join(["Matrix dim =", truncatedfockdim]));
     %
    % Show memory usage of the 5 largest variables.
    %
    bytecount = whos();
    varname = {bytecount.name};
    [bytecount, byteorder] = maxk([bytecount.bytes] / 1048576, 5);
    disp(join(["Memory used by variables [", varname{byteorder}, "] = [", bytecount, "] MBytes"]));

    disp('Start diagonalizing Hamiltonian matrix...');
    
    [nsigmaeigenstates, nsigmaspectra] = eig(ham, 'vector'); % NOTE : nsigmaeigenstates are column vectors
    % [nsigmaeigenstates, nsigmaspectra] = eigs(ham, truncatedfockdim);
    clear('ham');
    % Computing the matrix elements of the observables.
    if numel(varargin)
        if sun_symm
            testn(scidx) = sum(spinnocofig,'all');
        else
            testn{scidx} = spinnocofig;
        end
        if nn
            if ~sun_symm
                ninj_mean = cell(1,uint8(m*(m-1)/2)+m);
                pair_counter = 1;
		% different spin
                for idx = 1:m
                    for jdx = idx+1:m
                        if l>1
                            ninj_mean{pair_counter} = (basisvector(:,posi,idx) .* basisvector(:,posj,jdx))' * abs(nsigmaeigenstates).^2;
                        else
                            ninj_mean{pair_counter} = 0;
                        end
                        pair_counter = pair_counter + 1;
                    end
                end
		% same spin
		for idx=1:m
			if l>1
				ninj_mean{end + 1} = (basisvector(:,posi,idx).* basisvector(:,posj,idx))' * abs(nsigmaeigenstates).^2;
			else
				ninj_mean{end+1} = 0;
			end
			pair_counter = pair_counter + 1;
		end
                ninj{scidx} = ninj_mean;
                clear('ninj_mean');
            else
                cnnmatrix = 0;
                if l>1
                    for spinidx = 1:m
                        spina = repmat(spinidx, [1 m-1]);
                        spinb = [1:spinidx-1 spinidx+1:m];
                        cnnmatrix = cnnmatrix + sum(basisvector(:, posi, spina) .* (basisvector(:, posj, spina) - basisvector(:, posj, spinb)), 3);
                    end
                    ninjmatrix = permuteno(scidx) * (onsiteparticleno(:, posi) .* onsiteparticleno(:, posj))' * abs(nsigmaeigenstates).^2;
                    cnnmatrix = permuteno(scidx) * cnnmatrix' * abs(nsigmaeigenstates).^2;
                else
                    ninjmatrix = 0;
                    cnnmatrix = 0;
                end
                cnn{scidx} = cnnmatrix;
                ninj{scidx} = ninjmatrix;
                cnnmatrix = 0;
                ninjmatrix = 0;
            end   
        end
        if ni
            if ~sun_symm
                ni_mean = cell(1,m);
                for spin_idx = 1:m
                    ni_mean{spin_idx} = basisvector(:,:,spin_idx)' * abs(nsigmaeigenstates).^2;
                end
                nimatrix{scidx} = ni_mean;
            else
                nimatrix{scidx} = permuteno(scidx) * onsiteparticleno' * abs(nsigmaeigenstates).^2;
            end
        end
        
        if do
            if ~sun_symm
                do_mean = cell(1,uint8(m*(m-1)/2));
                for pair_idx = 1: uint8(m*(m-1)/2)
                    do_mean{pair_idx} = pair_count(pair_idx,:) * abs(nsigmaeigenstates).^2;
                end
                domatrix{scidx} = do_mean;
            else
                domatrix{scidx} = permuteno(scidx) * onsiteinteraction' * abs(nsigmaeigenstates).^2;
            end
        end
        if tri
            if ~sun_symm
                tri_mean = cell(1,uint8(m*(m-1)*(m-2)/6));
                for triplon_idx = 1: uint8(m*(m-1)*(m-2)/6)
                    tri_mean{triplon_idx} = triplon_count(triplon_idx,:) * abs(nsigmaeigenstates).^2;
                end
                trimatrix{scidx} = tri_mean;
            else
                trimatrix{scidx} = permuteno(scidx) * triplon_count' * abs(nsigmaeigenstates).^2;
            end
        end
    end
    spectra{scidx} = nsigmaspectra;
    
    % Clearing up memory.
    clear('nsigmaspectra','nsigmaeigenstates','tmatrix','basisvector','onsiteparticleno','onsiteinteraction');
end
%}
% Parallelize the large matrices
parfor scidx = 1:total_spin_permute_no
    spinnocofig = nsigma_sort{scidx};
    disp(join([scidx, "/", total_spin_permute_no, "nsigma species = [", spinnocofig, "] start diagonalizing"]));
    basisvector = cell(1, m); hmatrix = cell(1, m);
    dimlist = dimlist_arr_sort(scidx,:);
    fockdim = fockdimlist_sort(scidx);
    for spidx = 1:m
        basisvector{spidx} = subbasisvectorlist{spinnocofig(spidx) + 1};
        hmatrix{spidx} = hoppinglist{spinnocofig(spidx) + 1};
    end
     
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
    % clear('hmatrix', 'C');
    hmatrix = [];
    C = [];
    onsiteparticleno = sum(basisvector,3);
    %module for the interaction matrix and particle numbers
    if ~sun_symm
        onsiteinteraction = zeros(fockdim,1);
        pair_count = zeros(uint8(m*(m-1)/2),fockdim);
        triplon_count = zeros(uint8(m*(m-1)*(m-2)/6),fockdim);
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
                if numel(flavor_index)>2 % works for only 3 flavors. update for general n-flavors.
                    triplon_count(1,fockdim_idx) = triplon_count(1,fockdim_idx) + 1;
                end
            end
        end
        truncatedfockdim = fockdim;
        % Include truncation
        %
        %
    else
        onsiteinteraction = dot(onsiteparticleno, onsiteparticleno - 1, 2) / 2;
               
        % Basis-state truncation
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
        % p2_op = eq(onsiteparticleno,2);
        triplon_count = sum(onsiteparticleno .*(onsiteparticleno-1) .* (onsiteparticleno-2),2)/6;

    end
  
    ham = tmatrix;
    tmatrix = [];
    
    if ~truncatedfockdim
        disp('Hamiltionian matrix empty, skip to the next sector.');
        continue
    end

    % NOTE : the matrix is only upper triangular
    ham = ham + ham'; % ' is already Hermitian conjugate
    if sun_symm
        ham = ham + u*spdiags(onsiteinteraction, 0, truncatedfockdim, truncatedfockdim);
    else
        ham = ham + spdiags(onsiteinteraction, 0, truncatedfockdim, truncatedfockdim);
    end
    ham = full(ham);
 
    disp(join(["Matrix dim =", truncatedfockdim]));
     %
    % Show memory usage of the 5 largest variables.
    %
    % bytecount = whos();
    % varname = {bytecount.name};
    % [bytecount, byteorder] = maxk([bytecount.bytes] / 1048576, 5);
    % disp(join(["Memory used by variables [", varname{byteorder}, "] = [", bytecount, "] MBytes"]));

    disp('Start diagonalizing Hamiltonian matrix...');
    
    [nsigmaeigenstates, nsigmaspectra] = eig(ham, 'vector'); % NOTE : nsigmaeigenstates are column vectors
    % [nsigmaeigenstates, nsigmaspectra] = eigs(ham, truncatedfockdim);
    ham = []; % clearing memory

    % Computing the matrix elements of the observables.
    if numel(varargin)
        if sun_symm
            testn(scidx) = sum(spinnocofig,'all');
        else
            testn{scidx} = spinnocofig;
        end
        if nn
            if ~sun_symm
                ninj_mean = {};
                for idx = 1:m
                    for jdx = idx+1:m
                        if l>1
                            ninj_mean{end + 1} = (basisvector(:,posi,idx) .* basisvector(:,posj,jdx))' * abs(nsigmaeigenstates).^2;
                        else
                            ninj_mean{end + 1} = 0;
                        end
                    end
                end
                ninj{scidx} = ninj_mean;
                ninj_mean = {};
            else
                cnnmatrix = 0;
                if l>1
                    for spinidx = 1:m
                        spina = repmat(spinidx, [1 m-1]);
                        spinb = [1:spinidx-1 spinidx+1:m];
                        cnnmatrix = cnnmatrix + sum(basisvector(:, posi, spina) .* (basisvector(:, posj, spina) - basisvector(:, posj, spinb)), 3);
                    end
                    ninjmatrix = permuteno(scidx) * (onsiteparticleno(:, posi) .* onsiteparticleno(:, posj))' * abs(nsigmaeigenstates).^2;
                    cnnmatrix = permuteno(scidx) * cnnmatrix' * abs(nsigmaeigenstates).^2;
                else
                    ninjmatrix = 0;
                    cnnmatrix = 0;
                end
                cnn{scidx} = cnnmatrix;
                ninj{scidx} = ninjmatrix;
                cnnmatrix = 0;
                ninjmatrix = 0;
            end   
        end
        if ni
            if ~sun_symm
                ni_mean = cell(1,m);
                for spin_idx = 1:m
                    ni_mean{spin_idx} = basisvector(:,:,spin_idx)' * abs(nsigmaeigenstates).^2;
                end
                nimatrix{scidx} = ni_mean;
            else
                nimatrix{scidx} = permuteno(scidx) * onsiteparticleno' * abs(nsigmaeigenstates).^2;
            end
        end
        
        if do
            if ~sun_symm
                do_mean = cell(1,uint8(m*(m-1)/2));
                for pair_idx = 1: uint8(m*(m-1)/2)
                    do_mean{pair_idx} = pair_count(pair_idx,:) * abs(nsigmaeigenstates).^2;
                end
                domatrix{scidx} = do_mean;
            else
                domatrix{scidx} = permuteno(scidx) * onsiteinteraction' * abs(nsigmaeigenstates).^2;
            end
        end
        if tri
            if ~sun_symm
                tri_mean = cell(1,uint8(m*(m-1)*(m-2)/6));
                for triplon_idx = 1: uint8(m*(m-1)*(m-2)/6)
                    tri_mean{triplon_idx} = triplon_count(triplon_idx,:) * abs(nsigmaeigenstates).^2;
                end
                trimatrix{scidx} = tri_mean;
            else
                trimatrix{scidx} = permuteno(scidx) * triplon_count' * abs(nsigmaeigenstates).^2;
            end
        end
    end
    spectra{scidx} = nsigmaspectra;
    
    % Clearing up memory.
    nsigmaspectra = []; nsigmaeigenstates = [];
    tmatrix = []; basisvector = []; onsiteparticleno = []; 
    onsiteinteraction = [];
end

timer_count = toc(timer_count);

%saving the outputs in a struct
output = struct();
output.timer_count = timer_count;
output.spectra = spectra;
output.nsigmapermute = permuteno;

if numel(varargin)
    output.testn = testn;
    if ni
        output.nimatrix = nimatrix;
    end
    if ni2
        output.ni2matrix = ni2matrix;
    end
    if do
        output.domatrix = domatrix;
    end
    if tri
        output.trimatrix = trimatrix;
    end
    if nn
        output.ninj = ninj;
    end

save(filename, 'output', '-v7.3'); % save all result of the ED
disp(join([key, "l=", l, "t=", double(t), "u=", double(u), "n=", n, "m=", m, "D=", dno, "ED finish"]));
end

end


