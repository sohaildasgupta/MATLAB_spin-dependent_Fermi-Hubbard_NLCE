function accum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, varargin)
    %
    % LDA_measure - calculate thermodynamic observables from given energy spectrum and eigenstate expectation values.
    %
    % SYNTAX:
    %         accum = LDA_measure(spectra, deltamu, T, nsigmapermute, testn, matrix)
    %
    % INPUT ARGUMENTS:
    % spectra : energy spectra
    % deltamu : the orignal chemical potential subtract the chemical potential where you want to calculate observables
     if nargin == 6
        % Set default value if the argument is not provided
        matrix = varargin{1};
        sq = false;
    else
        % Use the provided value
        matrix = varargin{1};
        sq = true;
    end

    baslen = length(spectra);
    accum = 0;
    partitionfunction = 0;
    bias = min(cellfun(@(x) x(1), spectra) + cellfun(@(x) dot(x, deltamu), testn)); % Find GS energy and shift the energy zero to prevent partitionfunction blowup

    if iscell(matrix)
        for idx = 1:baslen
            nsigmaspectra = spectra{idx} + dot(deltamu , testn{idx});
            expspectra = exp(-(nsigmaspectra - bias) / T); 
            if ~sq
                accum = accum + matrix{idx} * expspectra;
            else
                accum = accum + matrix{idx}.^2 * expspectra;
            end
            partitionfunction = partitionfunction + nsigmapermute(idx) * sum(expspectra, 'all');
        end
        accum = accum / partitionfunction;
    else
        for idx = 1:baslen
            permuteno = nsigmapermute(idx);
            nsigmaspectra = spectra{idx} + dot(deltamu , testn{idx});
            expspectra = exp(-(nsigmaspectra - bias) / T);
            if strcmp(matrix,'energy_square')
                accum = accum + permuteno * (spectra{idx}.^2)' *expspectra;
            else
                accum = accum + permuteno * nsigmaspectra' * expspectra;
            end
            partitionfunction = partitionfunction + permuteno * sum(expspectra, 'all');
        end

        accum = accum / partitionfunction;

        if strcmp(matrix, 'entropy')% entropy, otherwise energy only
            accum = accum/T + log(partitionfunction) - bias/T;
        end
    end
end
