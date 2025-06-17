function threshold = entropythresh(densityMap)
    % Normalize the density map to the range [0, 1]
    densityMap = double(densityMap);
    densityMap = (densityMap - min(densityMap(:))) / (max(densityMap(:)) - min(densityMap(:)));
    
    % Compute the histogram of the normalized density map
    numBins = 256; % Number of bins for the histogram
    [counts, binEdges] = histcounts(densityMap, numBins);
    binLevels = binEdges(1:end-1) + diff(binEdges) / 2; % Get bin centers

    % Convert histogram counts to probabilities
    probabilities = counts / sum(counts);

    % Initialize maximum entropy variables
    maxEntropy = -inf;
    bestThreshold = 0;

    % Loop through all possible thresholds
    for i = 1:numBins
        % Background (below threshold)
        P_background = probabilities(1:i);
        weights_background = P_background / sum(P_background);
        entropyBackground = -sum(weights_background .* log2(weights_background + eps)); % Add eps to avoid log(0)

        % Foreground (above threshold)
        P_foreground = probabilities(i+1:end);
        weights_foreground = P_foreground / sum(P_foreground);
        entropyForeground = -sum(weights_foreground .* log2(weights_foreground + eps)); % Add eps to avoid log(0)

        % Total entropy
        totalEntropy = entropyBackground + entropyForeground;

        % Check if this threshold is better
        if totalEntropy > maxEntropy
            maxEntropy = totalEntropy;
            bestThreshold = binLevels(i);
        end
    end

    % Convert threshold back to the original scale of densityMap
    minDensity = min(densityMap(:));
    maxDensity = max(densityMap(:));
    threshold = bestThreshold * (maxDensity - minDensity) + minDensity;
end
