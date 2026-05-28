function output = ED_solver_parallel(graph, t, u, m, filename ,varargin)
%
% ED_solver - solve ED of a given graph or lattice size with different parameters and arguments
% Number conservation symmetry is used to block diagonalize the Hamiltonian matrix, and the function 
% can calculate the spectrum and observables for each block. The function can also automatically 
% generate the graph for a given lattice size with periodic boundary condition. 
% The function is parallelized for large matrices, and the output is saved in a struct for later use.
%
% SYNTAX :n
%         [spectra, nsigmapermute] = ED_solver(graph, t, u, m)
%         [spectra, nsigmapermute, testn, matrix] = ED_solver(graph, t, u, m, observables)
%
% INPUT ARGUMENTS :
% graph : input graph. If is a graph object, then the funciton reads he object; if is a lattice size vector [sz1, sz2, sz3], then the function automatically generates a graph of that dimension with periodci boundary condtion. NOTE: size vector dimension should equal to space dim, 1D 6x1 chain should bw wrtten as [6] instead of [6 1].
% t, u : parameters of FHM
% m : # of spin flavor species
%
% observables : 'density', 'doublon', 'triplon', 'nearest_pair'
% NOTE : 1. 'nearest_pair' contains two inputs & one output
%            inputs : list of postion i & postion j
%            output : state expectation values of ninj (density-density correlation between site i and j between same and different spin flavors)
%
% OUTPUT ARGUMENTS :
% spectra : energy spectrum of calculated system
% nsigmapermute : information of spin permutation
% testn : total ptcl # of the sector, stored for calculating results changing mu & T without diagonalizing Hamiltonian again i.e.: thermodynamics calculation
% matrix : state expectation values of every eigenstates, stored for calculating results changing mu & T without diagonalizing Hamiltonian again i.e.: thermodynamics calculation
%
% NOTE : 1. The basis states are ordered by spin flavor first, then by site index. e.g.: if |ud> = C^dagger_1u C^dagger_2d|0>, then |du> = - C^dagger_2u C^dagger_1d|0>, so the SU(2) singlet |ud>-|du>=C^dagger_1u C^dagger_2d|0> + C^dagger_2u C^dagger_1d|0>. As a result, the SU(2) spin singlet has the same phase factor for both components.
%
% Use graph key for NLCE summation.
% Authors : HT Wei and Sohail Dasgupta, Rice
%
    [ni, do, tri, nn] = deal(false);
    itemlist = {};

    timer_count = tic;

    l = numnodes(graph);

    if numel(varargin)
        input = varargin;
        for idx = 1:numel(input)
            arg = input{idx};
            if strcmp(arg, 'density')
                ni = true;
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
            elseif strcmp(arg, 'nearest_pair')
                nn = true;
                itemlist{end + 1} = arg;
                if l>1
                    posi = input{idx + 1};
                    posj = input{idx + 2};
                end
                continue
            end
        end
    end

    key = key_gen(graph);
    disp(join([key, "l=", l, "t=", double(t), "u=", double(u), "m=", m, "ED start"]));
    [side1, side2] = findedge(graph);
    orderlist = [side1 side2];

    disp(join(["Output options :", itemlist]));

    nsigmatuples = unique(sort(permn(0:l, m), 2), 'row'); 

    % Generate all the spin permutations
    nsigma = {}; 
    permuteno = [];

    for spinnocofig = nsigmatuples' %matlab iterates over columns
        spin_permutations = unique(perms(spinnocofig), 'row'); 
        for permute_idx = 1:size(spin_permutations,1)
            nsigma{end + 1} = spin_permutations(permute_idx,:);
            permuteno(end + 1) = 1;
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
    testn = cell(1,total_spin_permute_no);
    if ni
        nimatrix = cell(1,total_spin_permute_no);
    end
    if do
        domatrix = cell(1,total_spin_permute_no);
    end
    if tri
        trimatrix = cell(1,total_spin_permute_no);
    end
    if nn
        ninj = cell(1,total_spin_permute_no);
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

    % Parallelly diagonalize the Hamiltonian matrix for each sector, and compute the observables. 
    parfor scidx = 1:total_spin_permute_no
        spinnocofig = nsigma{scidx};
        disp(join([scidx, "/", total_spin_permute_no, "nsigma species = [", spinnocofig, "] start diagonalizing"]));
        basisvector = cell(1, m); hmatrix = cell(1, m);
        dimlist = dimlist_arr(scidx,:);
        fockdim = fockdimlist(scidx);
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
        % Clearing up memory for the hopping matrix.
        hmatrix = [];
        C = [];
        onsiteparticleno = sum(basisvector,3);
        
        % Module for the interaction matrix and particle numbers
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

        ham = tmatrix;
        tmatrix = [];
        
        % NOTE : the matrix is only upper triangular
        ham = ham + ham'; % ' is already Hermitian conjugate
        ham = ham + spdiags(onsiteinteraction, 0, fockdim, fockdim);
        ham = full(ham);
    
        disp(join(["Matrix dim =", fockdim]));
        disp('Start diagonalizing Hamiltonian matrix...');
        
        [nsigmaeigenstates, nsigmaspectra] = eig(ham, 'vector'); % NOTE : nsigmaeigenstates are column vectors
        ham = []; % clearing memory

        % Computing the matrix elements of the observables.
        if numel(varargin)
            testn{scidx} = spinnocofig;
            if nn
                ninj_mean = {};
                %different spin
                for idx = 1:m
                    for jdx = idx+1:m
                        if l>1
                            ninj_mean{end + 1} = (basisvector(:,posi,idx) .* basisvector(:,posj,jdx))' * abs(nsigmaeigenstates).^2;
                        else
                            ninj_mean{end + 1} = 0;
                        end
                    end
                end
                %same spin
                for idx=1:m
                    if l>1
                        ninj_mean{end + 1} = (basisvector(:,posi,idx) .* basisvector(:,posj,idx))' * abs(nsigmaeigenstates).^2;
                    else
                        ninj_mean{end + 1} = 0;
                    end
                end
                ninj{scidx} = ninj_mean;
                ninj_mean = {};
            end
            if ni
                ni_mean = cell(1,m);
                for spin_idx = 1:m
                    ni_mean{spin_idx} = basisvector(:,:,spin_idx)' * abs(nsigmaeigenstates).^2;
                end
                nimatrix{scidx} = ni_mean;
            end
            if do
                do_mean = cell(1,uint8(m*(m-1)/2));
                for pair_idx = 1: uint8(m*(m-1)/2)
                    do_mean{pair_idx} = pair_count(pair_idx,:) * abs(nsigmaeigenstates).^2;
                end
                domatrix{scidx} = do_mean;
            end
            if tri
                tri_mean = cell(1,uint8(m*(m-1)*(m-2)/6));
                for triplon_idx = 1: uint8(m*(m-1)*(m-2)/6)
                    tri_mean{triplon_idx} = triplon_count(triplon_idx,:) * abs(nsigmaeigenstates).^2;
                end
                trimatrix{scidx} = tri_mean;
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
        if do
            output.domatrix = domatrix;
        end
        if tri
            output.trimatrix = trimatrix;
        end
        if nn
            output.ninj = ninj;
        end

    save(filename, 'output'); % save all result of the ED
    disp(join([key, "l=", l, "t=", double(t), "u=", double(u), "m=", m, "ED finish"]));
    end
end