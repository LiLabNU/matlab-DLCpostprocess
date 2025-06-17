function H_grid = path_entropy_grid9(x, y, nbins, minX, maxX, minY, maxY)
    % Divide arena into 3x3 grid, compute entropy in each region
    %
    % Inputs:
    %   x, y    - position vectors
    %   nbins   - number of bins per sub-region (suggest small like 5)
    %   minX, maxX, minY, maxY - global limits
    %
    % Output:
    %   H_grid  - 3x3 matrix of entropy values (in bits)

    % Define grid divisions
    x_edges = linspace(minX, maxX, 4); % 3 regions → 4 edges
    y_edges = linspace(minY, maxY, 4);

    H_grid = NaN(3,3); % 3x3 output

    for i = 1:3
        for j = 1:3
            % Define subregion bounds
            x_min = x_edges(i); x_max = x_edges(i+1);
            y_min = y_edges(j); y_max = y_edges(j+1);

            % Get indices of positions within subregion
            idx = x >= x_min & x < x_max & y >= y_min & y < y_max;

            % Subsample trajectory in this subregion
            x_sub = x(idx);
            y_sub = y(idx);

            % Skip if no data
            if numel(x_sub) < 2
                H_grid(j,i) = NaN;  % (j,i) for rows=Y, cols=X
                continue;
            end

            % Compute entropy using same logic
            edgesX = linspace(x_min, x_max, nbins+1);
            edgesY = linspace(y_min, y_max, nbins+1);
            counts = histcounts2(x_sub, y_sub, edgesX, edgesY);

            p = counts(:);
            p = p / sum(p);
            p(p == 0) = [];

            H_grid(j,i) = -sum(p .* log2(p));
        end
    end
end
