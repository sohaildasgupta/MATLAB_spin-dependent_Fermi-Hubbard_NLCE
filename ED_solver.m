function output = ED_solver(graph, t, u, m, filename, varargin)
%
% ED_solver - solve ED of a given graph or lattice size with different parameters and arguments.
% Number conservation symmetry is used to block diagonalize the Hamiltonian matrix, and the function 
% can calculate the spectrum and observables for each block. The function can also automatically 
% generate the graph for a given lattice size with periodic boundary condition. 
% For parallelization, please use ED_solver_parallel.m instead, which has the same syntax and input/output arguments as this function.
%
% SYNTAX :n
%         [spectra, nsigmapermute] = ED_solver(graph, t, u, m)
%         [spectra, nsigmapermute, testn, matrix] = ED_solver(graph, t, u, m, observables)   
%          (saves the matrix elements to file.)
%
% INPUT ARGUMENTS :
% graph : input graph. If is a graph object, then the funciton reads he object; if is a lattice size vector [sz1, sz2, sz3], then the function automatically generates a graph of that dimension with periodci boundary condtion. NOTE: size vector dimension should equal to space dim, 1D 6x1 chain should bw wrtten as [6] instead of [6 1].
% t, u : parameters of FHM
% m : # of spin flavor species
%
% observables : 'density', 'doublons', 'triplons', 'nearest_pair'
% NOTE : 1. 'nearest_pair' contains two inputs & two outputs
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
% Authors : HT Wei and Sohail Dasgupta @ Rice
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
                break
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
    total_spin_permute_no = 0; %Counts the total number of spin permutations
    for spinnocofig = nsigmatuples' %matlab iterates over columns
        spin_permutations = unique(perms(spinnocofig), 'row'); 
        total_spin_permute_no = total_spin_permute_no + size(spin_permutations,1);
        nsigma{end + 1} = spin_permutations;
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
                if ridx < cidx % Only construct upper triangular part
                    connection = sort(swap);
                    % To give the right sign for correction indexing of basis states
                    v = -(-1)^sum(startstate((connection(1) + 1):(connection(end) - 1))) * t;
                    [idx(end + 1), jdx(end + 1), value(end + 1)] = deal(ridx, cidx, v);
                end
            end
        end
        hoppinglist{end + 1} = sparse(idx, jdx, value, dim, dim);
    end
    
    hoppinglist = {};
    for subbasisvector = subbasisvectorlist
        hfunction();
    end
    
    spectra = {}; nsigmapermute = [];     
    testn = {};
    if ni
        nimatrix = {};
    end
    if do
        domatrix = {};
    end
    if tri
        trimatrix = {};
    end
    if nn
        ninj = {};
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
        permuteno = 1; 
        spinconfig_len = size(spin_configs,1);

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
    
            %module for the interaction matrix and particle numbers
            onsiteparticleno = sum(basisvector,3);
           
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
                    if numel(flavor_index)>2 % works for only 3 flavors. 
                        triplon_count(1,fockdim_idx) = triplon_count(1,fockdim_idx) + 1;
                    end
                end
            end
           
           
            ham = tmatrix;
            clear('tmatrix');
            if ~fockdim
                disp('Hamiltionian matrix empty, skip to the next sector.');
                continue
            end

            % NOTE : the matrix is only upper triangular
            ham = ham + ham'; % ' is already Hermitian conjugate
            ham = ham + spdiags(onsiteinteraction, 0, fockdim, fockdim);
            ham = full(ham);
            
            disp(join(["Matrix dim =", fockdim]));
            varmem();

            disp('Start diagonalizing Hamiltonian matrix...');
            % tic
            
            [nsigmaeigenstates, nsigmaspectra] = eig(ham, 'vector'); % NOTE : nsigmaeigenstates are column vectors
            
            % toc
            clear('ham');

            if numel(varargin)
                testn{end + 1} = spinnocofig;
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
                    ninj{end + 1} = ninj_mean;
                    clear('ninj_mean');                       
                end
                if ni                        
                    ni_mean = {};
                    for spin_idx = 1:m
                        ni_mean{end+1} = basisvector(:,:,spin_idx)' * abs(nsigmaeigenstates).^2;
                    end
                    nimatrix{end + 1} = ni_mean;
                end
                if do
                    do_mean = {};
                    for pair_idx = 1: uint8(m*(m-1)/2)
                        do_mean{end+1} = pair_count(pair_idx,:) * abs(nsigmaeigenstates).^2;
                    end
                    domatrix{end + 1} = do_mean;
                end
                if tri
                    tri_mean = {};
                    for triplon_idx = 1: uint8(m*(m-1)*(m-2)/6)
                        tri_mean{end+1} = triplon_count(triplon_idx,:) * abs(nsigmaeigenstates).^2;
                    end
                    trimatrix{end + 1} = tri_mean;                        
                end
            end
            spectra{end + 1} = nsigmaspectra;
            nsigmapermute(end + 1) = permuteno;
            clear('nsigmaspectra','nsigmaeigenstates','tmatrix','basisvector',...
             'onsiteparticleno', 'originalonsitptclno', 'onsiteinteraction');            
        end   
    end
    
    timer_count = toc(timer_count);
    
    % saving the outputs in a struct
    output = struct();
    output.timer_count = timer_count;
    output.spectra = spectra;
    output.nsigmapermute = nsigmapermute;
    
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
    
    save(filename, 'output', '-v7'); % save all result of the ED
    
    disp(join([key, "l=", l, "t=", double(t), "u=", double(u), "m=", m, "ED finish"]));
    end
end