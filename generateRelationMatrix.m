function R = generateRelationMatrix(Y, task_num, target)
    % Generates a relation matrix based on the input parameters.
    % Inputs:
    %   Y - Input data vector
    %   task_num - Total number of tasks
    %   target - Number of targets
    % Output:
    %   R - Relation matrix

% Calculate the number of time segments
    TT = task_num / target;
% Initialize R for the time smoothness constraints
    R = zeros(TT, TT - 1);
    R(1:(TT + 1):end) = 1;
    R(2:(TT + 1):end) = -1;

 % Create the RR1T matrix
    RR1T = R;

    % Expand RR1T for the number of targets
    R = kron(eye(target), RR1T);

    % Initialize identity matrix M1
    M1 = eye(TT);

% Prepare correlation storage
    cor_target = cell(TT, 1);
    corr_target = zeros(target);
    
    for i = 1:TT
        u=1;
        for  j = 0 : target-1
            y_subset(u) = Y(j*TT+i);  
            u=u+1;
        end
        cor_target{i}=y_subset;
    end

   % Initialize matrix M
    M = zeros(target * TT);
% Populate the correlation matrix and M
    for i = 1:target
        for j = 1:target
            if i == j
                % If i equals j, set diagonal value
                corr_target(i, j) = 1 * target;
                M((i-1)*TT+1:i*TT, (j-1)*TT+1:j*TT) = M1 * corr_target(i, j);
            else
                % Calculate the correlation for off-diagonal elements
                corr_target_t = zeros(TT, 1); % Preallocate for performance
                for t = 1:TT
                    corr_target(i, j) = -corr(cor_target{t}{i}, cor_target{t}{j});
                    corr_target_t(t) = corr_target(i, j);
                end
                
                cor_diag = diag(corr_target_t); % Create a diagonal matrix
                M((i-1)*TT+1:i*TT, (j-1)*TT+1:j*TT) = M1 .* cor_diag;
            end
        end
    end

 % Final relation matrix calculation
    R = M * R;
    R = R';
end
