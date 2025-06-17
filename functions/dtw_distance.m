function distance = dtw_distance(traj1, traj2)
    % traj1 and traj2 are [N x 2] matrices of (x,y) points
    N1 = size(traj1, 1);
    N2 = size(traj2, 1);
    D = inf(N1+1, N2+1);
    D(1,1) = 0;
    for i = 1:N1
        for j = 1:N2
            cost = norm(traj1(i,:) - traj2(j,:)); % Euclidean distance
            D(i+1, j+1) = cost + min([D(i, j+1), D(i+1, j), D(i, j)]);
        end
    end
    distance = D(N1+1, N2+1);
end