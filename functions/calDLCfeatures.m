function results = calDLCfeatures(data, frameRate,shouldZscore)
    % Extract coordinates from the data structure
    nose_x = data.nose_1;
    nose_y = data.nose_2;
    head_x = data.head_1;
    head_y = data.head_2;
    body_x = data.body_1;
    body_y = data.body_2;

    % Coordinates for interaction points (specify or extract these based on your data)
    sucrosePort_x = data.sucroseport_1;
    sucrosePort_y = data.sucroseport_2;
    nosePortBot_x = data.botnoseport_1;
    nosePortBot_y = data.botnoseport_2;
    nosePortTop_x = data.topnoseport_1;  
    nosePortTop_y = data.topnoseport_2;

    % Calculate distances
    distanceToSucrosePort = sqrt((nose_x - sucrosePort_x).^2 + (nose_y - sucrosePort_y).^2);
    distanceToNosePortBot = sqrt((nose_x - nosePortBot_x).^2 + (nose_y - nosePortBot_y).^2);
    distanceToNosePortTop = sqrt((nose_x - nosePortTop_x).^2 + (nose_y - nosePortTop_y).^2);

    % Calculate head angles
    vector_head_to_nose = [nose_x - head_x, nose_y - head_y];
    vector_body_to_head = [head_x - body_x, head_y - body_y];
    dot_product = sum(vector_head_to_nose .* vector_body_to_head, 2);
    magnitude_head_to_nose = sqrt(sum(vector_head_to_nose.^2, 2));
    magnitude_body_to_head = sqrt(sum(vector_body_to_head.^2, 2));
    cos_theta = dot_product ./ (magnitude_head_to_nose .* magnitude_body_to_head);
    head_angles = acosd(cos_theta);  % angles in degrees

    % Calculate velocity and acceleration
    deltaTime = 1 / frameRate;  % Time between frames in seconds
    velocity_x = diff(nose_x) / deltaTime;
    velocity_y = diff(nose_y) / deltaTime;
    velocity = sqrt(velocity_x.^2 + velocity_y.^2);  % Magnitude of velocity
    acceleration_x = diff(velocity_x) / deltaTime;
    acceleration_y = diff(velocity_y) / deltaTime;
    acceleration = sqrt(acceleration_x.^2 + acceleration_y.^2);  % Magnitude of acceleration

    % Calculate turning angles
    angles = calculateTurningAngles(nose_x, nose_y);

    if shouldZscore
        distanceToSucrosePort = zscore(distanceToSucrosePort);
        distanceToNosePortBot = zscore(distanceToNosePortBot);
        distanceToNosePortTop = zscore(distanceToNosePortTop);
        head_angles = zscore(head_angles);
        velocity = zscore(velocity);
        acceleration = zscore(acceleration);
        angles = zscore(angles);
    end

    % Store results in a structure
    results.distanceToSucrosePort = distanceToSucrosePort;
    results.distanceToNosePortBot = distanceToNosePortBot;
    results.distanceToNosePortTop = distanceToNosePortTop;
    results.headAngles = head_angles;
    results.velocity = [velocity; velocity(end)];  % Append NaN to match the length of original data
    results.acceleration = [acceleration; acceleration(end); acceleration(end)];  % Append two NaNs
    results.turningAngles = [angles(1); angles; angles(end)];  % Append NaN to align lengths
end

function turningAngles = calculateTurningAngles(x, y)
    % Calculate turning angles from x and y coordinates
    x = x(:);
    y = y(:);
    dx = diff(x);
    dy = diff(y);
    vectorLengths = sqrt(dx.^2 + dy.^2);
    unitVectorsX = dx ./ vectorLengths;
    unitVectorsY = dy ./ vectorLengths;
    dotProducts = max(min(unitVectorsX(1:end-1) .* unitVectorsX(2:end) + ...
                  unitVectorsY(1:end-1) .* unitVectorsY(2:end), 1), -1);
    anglesRadians = acos(dotProducts);
    turningAngles = rad2deg(anglesRadians);
end
