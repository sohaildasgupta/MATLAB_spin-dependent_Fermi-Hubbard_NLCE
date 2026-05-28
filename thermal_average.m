function accum = thermal_average(spectra, deltamu, T, nsigmapermute, testn, varargin)
    %
    % thermal_average - calculate thermodynamic observables from given energy spectrum and eigenstate expectation values.
    %
    % SYNTAX:
    %         accum = thermal_average(spectra, deltamu, T, nsigmapermute, testn, matrix)
    %
    % INPUT ARGUMENTS:
    % spectra : energy spectra
    % deltamu : the orignal chemical potential subtract the chemical potential where you want to calculate observables
    % T : temperature
    % nsigmapermute : the number of states in each symmetry sector
    % testn : the quantum numbers of each symmetry sector
    % matrix : the observable you want to calculate, can be 'energy','entropy', or the matrix elements of the observable in each symmetry sector. 
    matrix = varargin{1};
    baslen = length(spectra);
    accum = 0;
    partitionfunction = 0;
    
    % Find GS energy and shift the energy zero to prevent partitionfunction blowup
    bias = min(cellfun(@(x) x(1), spectra) + cellfun(@(x) dot(x, deltamu), testn)); 

    if iscell(matrix)
        for idx = 1:baslen
            nsigmaspectra = spectra{idx} + dot(deltamu , testn{idx});
            expspectra = exp(-(nsigmaspectra - bias) / T); 
            accum = accum + matrix{idx} * expspectra;
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

        if strcmp(matrix, 'entropy')    % entropy, otherwise energy only
            accum = accum/T + log(partitionfunction) - bias/T;
        end
    end
end
