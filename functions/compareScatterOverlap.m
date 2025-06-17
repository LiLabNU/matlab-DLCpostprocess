function [overlapArea, jaccardIndex] = compareScatterOverlap(data1, data2, numGridPoints, showPlots, binaryOption)
    % compareScatterOverlap - Compares overlap between two scatter plots using KDE.
    % This version includes an option to convert density maps to binary and then compare overlaps.
    %
    % Syntax: [overlapArea, jaccardIndex] = compareScatterOverlap(data1, data2, numGridPoints, showPlots, binaryOption)
    %
    % Inputs:
    %   data1 - Nx2 matrix of xy coordinates for dataset 1
    %   data2 - Mx2 matrix of xy coordinates for dataset 2
    %   numGridPoints - Number of points in each dimension for the grid
    %   showPlots - Boolean flag to show plots or not
    %   binaryOption - Boolean flag to convert density maps to binary before overlap comparison
    %
    % Outputs:
    %   overlapArea - Scalar quantifying the estimated overlap area
    %   jaccardIndex - Scalar quantifying the Jaccard index of the overlap

    % Define grid for KDE
    minX = min([data1(:,1); data2(:,1)]);
    maxX = max([data1(:,1); data2(:,1)]);
    minY = min([data1(:,2); data2(:,2)]);
    maxY = max([data1(:,2); data2(:,2)]);
    x = linspace(minX, maxX, numGridPoints);
    y = linspace(minY, maxY, numGridPoints);
    [X, Y] = meshgrid(x, y);

    % Compute KDE for both datasets
    [F1, ~] = ksdensity(data1, [X(:), Y(:)]);
    [F2, ~] = ksdensity(data2, [X(:), Y(:)]);
    F1 = reshape(F1, size(X));
    F2 = reshape(F2, size(X));

    % Calculate binary overlap for Jaccard index
    if binaryOption
        B1 = F1 > 0;
        B2 = F2 > 0;
    else
        B1 = F1 > mean(F1(:)); % Use a threshold of the mean density
        B2 = F2 > mean(F2(:)); % Use a threshold of the mean density
    end



    overlapArea = sum(sum(min(B1, B2))); % Area of overlap in binary terms
    intersection = sum(sum(B1 & B2));
    union = sum(sum(B1 | B2));
    jaccardIndex = intersection / union;

    % Optional plot to visualize the overlap
    if showPlots
        figure;
        colormap('gray'); % Set colormap to grayscale

        subplot(1,3,1);
        surf(X, Y, double(B1), 'EdgeColor', 'none'); % Convert logical to double
        title('Binary Density of Dataset 1');
        view(2);
        colorbar;

        subplot(1,3,2);
        surf(X, Y, double(B2), 'EdgeColor', 'none'); % Convert logical to double
        title('Binary Density of Dataset 2');
        view(2);
        colorbar;

        subplot(1,3,3);
        contourf(X, Y, double(B1 & B2), 20, 'LineColor', 'none');
        title('Overlap of Binary Densities');
        colorbar;
    end

    % Return results
end
