function [xc, yc, R, resid] = fitCircle(x, y)
    % Number of data points
    n = length(x);

    % Formulate the system of equations
    A = [-2*x, -2*y, ones(n, 1)];
    b = -(x.^2 + y.^2);

    % Solve the system using least squares
    params = A\b;

    % Extract the circle parameters
    xc = params(1);
    yc = params(2);
    R = sqrt(params(3) + xc^2 + yc^2);

    % Calculate residuals
    fitted_circle = @(x, y) sqrt((x - xc).^2 + (y - yc).^2) - R;
    resid = sum(fitted_circle(x, y).^2);

    % Plot the fitted circle
    theta = linspace(0, 2*pi, 100);
    x_fit = xc + R*cos(theta);
    y_fit = yc + R*sin(theta);

    figure;
    scatter(x, y, 'b');
    hold on;
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    title('Fitted Circle');
    xlabel('X');
    ylabel('Y');
    legend('Data Points', 'Fitted Circle');
    axis equal;
    hold off;
end




