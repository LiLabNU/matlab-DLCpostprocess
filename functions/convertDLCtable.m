function bodyPartStruct = convertDLCtable(DLCdata, accuracyThreshold)
    % Extract all body part names dynamically
    allHeaders = DLCdata.Properties.VariableNames;
    bodyParts = unique(erase(allHeaders, ["_1", "_2"])); % Extract base body part names
    bodyParts(contains(bodyParts, "bodyparts", 'IgnoreCase', true)) = []; 
    
    % Initialize output structure
    bodyPartStruct = struct();
    
    % Process each body part
    for i = 1:length(bodyParts)
        bodyPart = bodyParts{i};
        x_col = bodyPart;        % X coordinate
        y_col = bodyPart + "_1"; % Y coordinate
        acc_col = bodyPart + "_2"; % Accuracy
        
        % Extract data
        x = DLCdata.(x_col);
        y = DLCdata.(y_col);
        accuracy = DLCdata.(acc_col);
        
        if ~strcmp(bodyPart,'sucroseport')
            % Find frames with accuracy below threshold
            badFrames = accuracy < accuracyThreshold;

            % Interpolate bad frames
            x(badFrames) = interp1(find(~badFrames), x(~badFrames), find(badFrames), 'linear', 'extrap');
            y(badFrames) = interp1(find(~badFrames), y(~badFrames), find(badFrames), 'linear', 'extrap');
        end

        % Store in structure (107980 x 2 array)
        bodyPartStruct.(bodyPart) = [x, y];
    end
end