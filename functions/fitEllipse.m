function [a, b, xc, yc, theta, resid] = fitEllipse(x, y)
    % Normalize data
    x = x(:);
    y = y(:);
    X = [x.^2, x.*y, y.^2, x, y, ones(size(x))];

    % Solve the system of equations
    [~, ~, V] = svd(X, 0);
    params = V(:, end);

    % Extract the ellipse parameters
    a = params(1);
    b = params(2);
    c = params(3);
    d = params(4);
    e = params(5);
    f = params(6);

    % Calculate the center of the ellipse
    xc = (b*e - 2*c*d) / (4*a*c - b^2);
    yc = (b*d - 2*a*e) / (4*a*c - b^2);

    % Calculate the angle of rotation of the ellipse
    theta = 0.5 * atan(b / (a - c));

    % Calculate the semi-major and semi-minor axes
    num = 2 * (a*e^2 + c*d^2 + f*b^2 - 2*b*d*e - a*c*f);
    den1 = (b^2 - a*c) * ((c - a) * sqrt(1 + 4*b^2 / ((a - c)^2)) - (c + a));
    den2 = (b^2 - a*c) * ((a - c) * sqrt(1 + 4*b^2 / ((a - c)^2)) - (c + a));
    a = sqrt(num / den1);
    b = sqrt(num / den2);

    % Calculate residuals
    fitted_ellipse = @(x, y) a*x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f;
    resid = sum(fitted_ellipse(x, y).^2);

    % Plot the fitted ellipse
    t = linspace(0, 2*pi, 100);
    X_fit = a * cos(t);
    Y_fit = b * sin(t);

    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    ellipse = R * [X_fit; Y_fit];

    figure;
    scatter(x, y, 'b');
    hold on;
    plot(ellipse(1, :) + xc, ellipse(2, :) + yc, 'r-', 'LineWidth', 2);
    title('Fitted Ellipse');
    xlabel('X');
    ylabel('Y');
    legend('Data Points', 'Fitted Ellipse');
    axis equal;
    hold off;
end


