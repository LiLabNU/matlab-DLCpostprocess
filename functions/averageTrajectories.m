function [avgTrajectoryX, avgTrajectoryY] = averageTrajectories(trajectories)
    % Assume trajectories is a cell array where each cell contains an [N x 2] array of (x, y) coordinates

    % Identify empty cells
    isEmpty = cellfun(@isempty, trajectories);

    % Remove empty cells
    trajectories = trajectories(~isEmpty);

    % Determine the common length - using the median length of trajectories
    lengths = cellfun(@(c) size(c, 1), trajectories);
    commonLength = round(median(lengths));

    % Initialize matrices to store interpolated trajectories
    interpolatedX = zeros(commonLength, numel(trajectories));
    interpolatedY = zeros(commonLength, numel(trajectories));

    % Interpolate each trajectory
    for i = 1:numel(trajectories)
        currentTrajectory = trajectories{i};
        oldPoints = linspace(0, 1, size(currentTrajectory, 1));
        newPoints = linspace(0, 1, commonLength);
        interpolatedX(:, i) = interp1(oldPoints, currentTrajectory(:, 1), newPoints, 'linear');
        interpolatedY(:, i) = interp1(oldPoints, currentTrajectory(:, 2), newPoints, 'linear');
    end

    % Calculate the average trajectory
    avgTrajectoryX = mean(interpolatedX, 2);
    avgTrajectoryY = mean(interpolatedY, 2);
end
