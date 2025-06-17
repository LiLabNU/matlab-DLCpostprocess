%% Update Sucrose Port Coordinates in DLC h5 Files
serverPath = "\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\PBS\LiPatel_Labs\Personal_Folders";
motherDir = "Talia\Behavior\MedPC_Data\Pavlov_pilot3\Day10";
folder = fullfile(serverPath, motherDir);

video_files = dir(fullfile(folder, '*.mp4')); % Update extension if needed
video_files = video_files(~contains({video_files.name}, 'DLC'));

for i = 1:length(video_files)
    [~, name, ~] = fileparts(video_files(i).name);
    video_path = fullfile(folder, video_files(i).name);

    h5_name = [name, 'DLC_resnet50_PavlovianOct11shuffle1_300000_filtered.h5'];
    h5_path = fullfile(folder, h5_name);

    if ~isfile(h5_path)
        warning('Corresponding h5 file not found for %s. Skipping...', video_files(i).name);
        continue;
    end

    % Read first frame
    v = VideoReader(video_path);
    firstFrame = readFrame(v);

    data = h5read(h5_path, '/df_with_missing/table');
    values = data.values_block_0';  % shape: [frames x features]

    % Read and parse column labels
    info = h5info(h5_path, '/df_with_missing');
    labels = info.Attributes(strcmp({info.Attributes.Name}, 'non_index_axes')).Value;

    % Extract readable column names
    tokens = regexp(labels, '\(V([^\n]+)\np\d+\nV([^\n]+)\np\d+\nV([^\n]+)\np\d+\)', 'tokens');
    colnames = cellfun(@(x) sprintf('%s_%s_%s', x{1}, x{2}, x{3}), tokens, 'UniformOutput', false);

    % Check label count matches data columns
    assert(length(colnames) == size(values, 2), 'Mismatch between column labels and data.');

    % Convert to table
    DLC_table = array2table(values, 'VariableNames', colnames);

    % Display first-frame x/y coordinates only
    xy_cols = contains(colnames, '_x') | contains(colnames, '_y');
    first_frame_xy = DLC_table(1, xy_cols);
    disp(first_frame_xy)

    % Extract head coordinates
    head_x_idx = find(colnames=="DLC_resnet50_PavlovianOct11shuffle1_300000_head_x");
    head_y_idx = find(colnames=="DLC_resnet50_PavlovianOct11shuffle1_300000_head_y");

    head_x = values(1, head_x_idx);
    head_y = values(1, head_y_idx);

    % Display frame with head coordinates
    figure; imshow(firstFrame);
    hold on;
    plot(head_x, head_y, 'ro', 'MarkerSize', 12, 'LineWidth', 2);
    text(head_x + 10, head_y, sprintf('(%.1f, %.1f)', head_x, head_y), 'Color','red', 'FontSize',12);
    title(sprintf('Click on Sucrose Port for %s', name), 'Interpreter','none');

    [x, y] = ginput(1);
    close;

    % Identify sucrose port column indices
    sucroseport_x_idx = find(colnames=="DLC_resnet50_PavlovianOct11shuffle1_300000_sucroseport_x");
    sucroseport_y_idx = find(colnames=="DLC_resnet50_PavlovianOct11shuffle1_300000_sucroseport_y");
    sucroseport_likelihood_idx = find(colnames=="DLC_resnet50_PavlovianOct11shuffle1_300000_sucroseport_likelihood");

    % Update coordinates for all frames
    values(:, sucroseport_x_idx) = x;
    values(:, sucroseport_y_idx) = y;
    if ~isempty(sucroseport_likelihood_idx)
        values(:, sucroseport_likelihood_idx) = 1;
    end

    % Save updated data
    h5write(h5_path, '/df_with_missing/table/values_block_0', values');

    fprintf('Updated sucrose port in %s\n', h5_name);
end

disp('All files processed.');