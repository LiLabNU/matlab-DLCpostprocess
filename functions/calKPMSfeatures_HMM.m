function results = calKPMSfeatures_HMM(data, frameRate)
   currAngle = rad2deg(data.heading);


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


