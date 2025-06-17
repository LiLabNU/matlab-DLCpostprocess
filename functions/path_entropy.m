function H = path_entropy(x, y, nbins, minX, maxX, minY, maxY)
    % PATH_ENTROPY computes the spatial entropy of a trajectory
    % 
    % Inputs:
    %   x, y      - Position vectors for one trial
    %   nbins     - Number of bins per dimension (scalar)
    %   minX/maxX/minY/maxY - Global spatial limits across all trials
    %
    % Output:
    %   H         - Shannon entropy (in bits)

    % Define bin edges
    edgesX = linspace(minX, maxX, nbins+1);
    edgesY = linspace(minY, maxY, nbins+1);

    % Compute 2D histogram
    counts = histcounts2(x, y, edgesX, edgesY);

    % Flatten to 1D and normalize to probabilities
    p = counts(:);
    p = p / sum(p);

    % Remove zeros to avoid log(0)
    p(p == 0) = [];

    % Compute entropy
    H = -sum(p .* log2(p));
end
