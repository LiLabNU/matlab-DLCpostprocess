function [results, Zresults] = calDLCfeatures_HMM(data, framerate, windowSize)
    % Extract coordinates from the data structure
    nose_x = data.nose_1;
    nose_y = data.nose_2;
    body_x = data.body_1;
    body_y = data.body_2;
    mouseCentroid = [body_x, body_y];

    % Coordinates for interaction points (specify or extract these based on your data)
    sucrosePort_x = data.sucroseport_1;
    sucrosePort_y = data.sucroseport_2;
    nosePortBot_x = data.botnoseport_1;
    nosePortBot_y = data.botnoseport_2;
    nosePortTop_x = data.topnoseport_1;
    nosePortTop_y = data.topnoseport_2;

    % Calculate distances
    distanceToSucrosePort = calculateEuclideanDistance([nose_x, nose_y], [sucrosePort_x(1), sucrosePort_y(1)]);
    distanceToNosePortBot = calculateEuclideanDistance([nose_x, nose_y], [nosePortBot_x(1), nosePortBot_y(1)]);
    distanceToNosePortTop = calculateEuclideanDistance([nose_x, nose_y], [nosePortTop_x(1), nosePortTop_y(1)]);

    % Remove outliers from distances
    distanceToSucrosePort = removeOutliers(distanceToSucrosePort);
    distanceToNosePortBot = removeOutliers(distanceToNosePortBot);
    distanceToNosePortTop = removeOutliers(distanceToNosePortTop);

    % Apply median filter to distances
    distanceToSucrosePort = medfilt1(distanceToSucrosePort, windowSize);
    distanceToNosePortBot = medfilt1(distanceToNosePortBot, windowSize);
    distanceToNosePortTop = medfilt1(distanceToNosePortTop, windowSize);

    % Calculate velocity
    velocity = zeros(length(body_x), 1);
    for i = 2:length(body_x)
        velocity(i) = norm([mouseCentroid(i, :)] - [mouseCentroid(i-1, :)]);
    end

    % Remove outliers from velocity
    velocity = removeOutliers(velocity);

    % Apply median filter to velocity
    velocity = medfilt1(velocity, windowSize);

    % Calculate acceleration
    dt = 1 / framerate;
    acceleration = gradient(velocity) ./ dt;

    % Remove outliers from acceleration
    acceleration = removeOutliers(acceleration);

    % Apply median filter to acceleration
    acceleration = medfilt1(acceleration, windowSize);

    % Store results in a structure
    results.distanceToSucrosePort = distanceToSucrosePort;
    results.distanceToNosePortBot = distanceToNosePortBot;
    results.distanceToNosePortTop = distanceToNosePortTop;
    results.velocity = velocity;
    results.acceleration = acceleration;

    % Normalize (z-score) the results
    distanceToSucrosePort = zscore(distanceToSucrosePort, 1);
    distanceToNosePortBot = zscore(distanceToNosePortBot, 1);
    distanceToNosePortTop = zscore(distanceToNosePortTop, 1);
    velocity = zscore(velocity, 1);
    acceleration = zscore(acceleration, 1);

    % Store z-scored results in a structure
    Zresults.distanceToSucrosePort = distanceToSucrosePort;
    Zresults.distanceToNosePortBot = distanceToNosePortBot;
    Zresults.distanceToNosePortTop = distanceToNosePortTop;
    Zresults.velocity = velocity;
    Zresults.acceleration = acceleration;
end

function distances = calculateEuclideanDistance(x, y)
    distances = zeros(size(x, 1), 1);
    for i = 1:size(x, 1)
        X = [y; x(i, :)];
        distances(i) = pdist(single(X), 'euclidean');
    end
end

