function [F, xgrid, ygrid, P] = estimateDensity(xyData, gridSize, bandwidth, plotDensity)
    % estimateDensity - Estimates the density of xy coordinates using KDE.
    % It also calculates and returns the probability distribution over the grid.
    %
    % Syntax: [F, xgrid, ygrid, P] = estimateDensity(xyData, gridSize, bandwidth, plotDensity)
    %
    % Inputs:
    %   xyData - Nx2 matrix of xy coordinates
    %   gridSize - Size of the grid (e.g., 100 for 100x100 grid points)
    %   bandwidth - 1x2 vector specifying the bandwidth for KDE
    %   plotDensity - Boolean, set to true to plot the density
    %
    % Outputs:
    %   F - Absolute density values on the grid
    %   xgrid - Vector of x grid points
    %   ygrid - Vector of y grid points
    %   P - Normalized probability distribution corresponding to F

    % Extract x and y data
    x = xyData(:,1);
    y = xyData(:,2);

    % Define grid over which to evaluate the density
    xgrid = linspace(min(x), max(x), gridSize);
    ygrid = linspace(min(y), max(y), gridSize);
    [Xgrid, Ygrid] = meshgrid(xgrid, ygrid);

    % Perform kernel density estimation
    [f, xi] = ksdensity([x, y], [Xgrid(:), Ygrid(:)], 'Bandwidth', bandwidth);

    % Reshape the output density into a 2D grid
    F = reshape(f, size(Xgrid));

    % Calculate the probability distribution by normalizing F
    totalArea = trapz(ygrid, trapz(xgrid, F, 2));  % Integral over the grid to find total area under F
    P = F / totalArea;  % Normalize F to get the probability distribution P

    % Plot the results if requested
    if plotDensity
        figure;
        subplot(1, 2, 1);
        surf(xgrid, ygrid, F, 'EdgeColor', 'none');
        title('Absolute Density Estimate');
        xlabel('X coordinate');
        ylabel('Y coordinate');
        zlabel('Density');
        colorbar;

        subplot(1, 2, 2);
        surf(xgrid, ygrid, P, 'EdgeColor', 'none');
        title('Probability Distribution');
        xlabel('X coordinate');
        ylabel('Y coordinate');
        zlabel('Probability');
        colorbar;
    end
end
