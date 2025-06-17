function [p, R2] = fitLinear(x, y, figNum)
% linearRegressionWithGoodnessOfFit performs linear regression on x and y data
% and plots the scatter plot along with the regression line.
%
% Inputs:
%   x - Vector of x data
%   y - Vector of y data
%   figNum - Figure number to plot the data (optional)
%
% Outputs:
%   p - Polynomial coefficients of the linear fit (degree 1)
%   R2 - Coefficient of determination (goodness of fit)

% Check if figNum is provided
if nargin < 3
    figNum = 0;
else
    figure(figNum);
end

% Perform linear regression
p = polyfit(x, y, 1); % Linear fit (degree 1)
yfit = polyval(p, x);

% Calculate R² (coefficient of determination)
ymean = mean(y);
SST = sum((y - ymean).^2);
SSR = sum((y - yfit).^2);
R2 = 1 - SSR/SST;

if figNum ~=0
    % Plot the original scatter plot
    scatter(x, y, 100, 'filled');
    hold on;

    % Plot the regression line
    plot(x, yfit, 'r-', 'LineWidth', 2);

    % Add labels and title with R² value
    xlabel('X-axis');
    ylabel('Y-axis');
    title(['Linear Regression (R^2 = ', num2str(R2, '%.2f'), ')']);
    legend('Data Points', 'Linear Fit');

    % Optionally, add R² as text annotation
    text(min(x), max(y), ['R^2 = ', num2str(R2, '%.2f')], 'FontSize', 12, 'Color', 'red');

    hold off;
end
end
