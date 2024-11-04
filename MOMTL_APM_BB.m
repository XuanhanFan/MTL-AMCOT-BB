%% FUNCTION MTL-AMCOT-BB
% Adaptive Multi-Cognitive Objective Temporal Task Approach for Predicting AD Progression" has been accepted as a regular paper at IEEE International Conference on Bioinformatics and Biomedicine 2024 (IEEE BIBM 2024).
%     APM-based Algorithm with Barzilai-Borwein Step Size
%% OBJECTIVE
%   argmin_W { sum_i^t (0.5 * norm (Y{i} - X{i}' * W(:, i))^2)
%            + rho1 * \|W\|_1 + rho2 * \|W*R\|_1  +
%            + rho3 * \|W\|_{2,1} }
%   where
%   rho1: sparse.
%   rho2: fused Lasso.
%   rho3: L2,1-norm.
%  R-encoding multiple cognitive goal scores and the relationship between time and task structure
%
%% INPUT
%   X: {d * n} * t - input matrix
%   Y: {n * 1} * t - output matrix
%   rho1: sparse.
%   rho2: fused Lasso.
%   rho3: L2,1-norm.
%
%% OUTPUT
%   W: model: d * t
%   funcVal: function value vector.
%
%% LICENSE
% This program is free software: developed based on Jiayu Zhou and Jieping Ye's MALSAR-master program. GNU General Public License
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. See <http://www.gnu.org/licenses/>.
%
% This project
%
% Last modified on 5 November 2024.
%
%% RELATED PAPERS
%
%   Adaptive Multi-Cognitive Objective Temporal Task Approach for Predicting AD Progression" has been accepted as a regular paper at IEEE International Conference on Bioinformatics and Biomedicine 2024 (IEEE BIBM 2024).
%     APM-based Algorithm with Barzilai-Borwein Step Size
%
%% RELATED FUNCTIONS
%   APM-based Algorithm, init_opts,generateRelationMatrix

function [W, funcVal,out] =  MTMTL_FISTABB(X, Y, target,rho1, rho2, rho3, opts)

if nargin < 6
    opts = [];
end

% Perform multi-dimensional transpose on X
X = multi_transpose(X);

% Initialize options. Set the opts variable to default options for use in the function.
opts = init_opts(opts); 
% % opts.ls: Indicates whether line search is used
% % opts.bb: Indicates whether BB step size is adopted 
% if ~isfield(opts, 'ls'); opts.ls = 1; end
% if ~isfield(opts, 'bb'); opts.bb = 1; end
opts.bb = 1;

task_num = length(X);
dimension = size(X{1}, 1);
funcVal = [];



R = generateRelationMatrix(Y, task_num, target);
% initial W
%W0 = zeros(dimension, task_num);
if opts.init==2
    W0 = zeros(dimension, task_num);
elseif opts.init== 0
    W0 = randn(dimension, task_num);
else

    if isfield(opts,'W0')
        W0=opts.W0;
        if (nnz(size(W0)-[dimension, task_num]))
            error('\n Check the input .W0');
        end
    else
        W0 = zeros(dimension, task_num);
    end
end

bFlag=0; % this flag tests whether the gradient step only changes a little


Wz = W0;
Wz_old = W0;
gWs_old = gradVal_eval(Wz_old); % Initialize gradient
t = 1;
t_old = 0;
tt = tic;
iter = 1;
gamma = 2; 
gamma_inc = 1/norm(X{1}, 'fro')^2;

bbover = 0; 
% Cval = 0; Q = 1; beta = 0.85; rhols = 1e-6;

while iter < opts.maxIter
    alpha = (t_old - 1) / t;

    Ws = (1 + alpha) * Wz - alpha * Wz_old;

    % Compute function value and gradients of the search point
    gWs = gradVal_eval(Ws);
    Fs = funVal_eval(Ws); % Objective function, cumulative error of least squares

    if iter > 1 && opts.bb
        % Calculate weight and gradient differences
        s = Wz - Wz_old;  % Weight difference
        y = gradVal_eval(Wz) - gWs_old;  % Gradient difference
        sty = abs(s(:)' * y(:));  % Step size calculation

        if sty > 0
            % Calculate step size using Barzilai-Borwein method
            if mod(iter, 2) == 0 
                gamma = (s(:)' * s(:)) / sty;  % Even iteration
            else 
                gamma = sty / (y(:)' * y(:));  % Odd iteration
            end
            % Optionally clamp gamma to a reasonable range
            % gamma = min(max(gamma, 0.5), 1e-5);
        end
    else
        % Adjust gamma if not using BB method
        gamma = gamma * 2;  % or gamma * 0.5; depending on your needs
    end
    
    % Perform line search to update weights
    [Wzp, Ffp, gamma] = lineSearch(Ws, gWs, gamma);
    Wz_old = Wz;  % Update old weights
    Wz = Wzp;  % Assign new weights
    gWs_old = gradVal_eval(Wz_old);  % Update old gradient
    
    % Record the function value for the current iteration
    funcVal = cat(1, funcVal, Ffp + nonsmooth_eval(Wz, rho1, rho2, rho3));


  
    
    if (bFlag)
        % fprintf('\n The program terminates as the gradient step changes the solution very small.');
        break;
    end
    
    % test stop condition.
    switch(opts.tFlag)
        case 0
            if iter>=2
                if (abs( funcVal(end) - funcVal(end-1) ) <= opts.tol)
                    break;
                end
            end
        case 1
            if iter>=2
                if (abs( funcVal(end) - funcVal(end-1) ) <=...
                        opts.tol* funcVal(end-1))
                    break;
                end
            end
        case 2
            if ( funcVal(end)<= opts.tol)
                break;
            end
        case 3
            if iter>=opts.maxIter
                break;
            end
    end
    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
end

W = Wzp;
out.tt = toc(tt); 
% private functions


function [Wp, Ff, gamma] = lineSearch(W, grad, gamma)
    % Initial Armijo parameter
    tau = 0.5;
    % c = 1e-4;  % Uncomment if using a constant for Armijo condition
    Wf = funVal_eval(W);  % Evaluate the function value at W

    while true
        % Project the updated W using the FGLasso projection
        Wp = FGLasso_projection(W - gamma * grad, rho1 * gamma, rho2 * gamma, rho3 * gamma);
        
        % Evaluate the function value at the projected W
        Ff = funVal_eval(Wp);
        delta_Wzp = Wp - W;  % Compute the change in W
        nrm_delta_Wzp = norm(delta_Wzp, 'fro')^2;  % Frobenius norm of the change
        r_sum = nrm_delta_Wzp;
        
        % Calculate the predicted function value based on the gradient
        Ffp_gamma = Wf + sum(sum(delta_Wzp .* grad)) + 1/gamma/2 * nrm_delta_Wzp;

        % Check for convergence conditions
        if (r_sum <= 1e-20)
            bFlag = 1;  % Indicates that the gradient step makes little improvement
            break;
        end
        
        if Ff <= Ffp_gamma
            break;  % Condition satisfied, exit the loop
        else
            gamma = tau * gamma;  % Reduce gamma if condition not satisfied
        end
        
        % Decrease alpha if needed
    end
end


    

    function [Wp] = FGLasso_projection (W, lambda_1, lambda_2, lambda_3 )
        % solve it in row wise (L_{2,1} is row coupled).
        % for each row we need to solve the proximal opterator
        % argmin_w { 0.5 \|w - v\|_2^2
        %            + lambda_1 * \|w\|_1 + lambda_2 * \|R * w\|_1  +
        %            + lambda_3 * \|w\|_2 }
        % NOTE: Here the R is t-1 * t, the outside R is t * t-1, and
        
        % Initialize Wp with the same size as W
        Wp = zeros(size(W));
        
        % Iterate over each row of W
        for i = 1:size(W, 1)
            v = W(i, :); % Extract the i-th row of W
            % Solve the proximal operator for the i-th row
            Wp(i, :) = FGLasso_projection_rowise(v, lambda_1, lambda_2, lambda_3)';
        end
end

% smooth part gradient.
    function [grad_W] = gradVal_eval(W)
    % Initialize grad_W based on the number of tasks
    grad_W = zeros(size(W));
    
    % Use parallel computation if opts.pFlag is true
    if opts.pFlag
        parfor i = 1:task_num
            grad_W(:, i) = X{i} * (X{i}' * W(:, i) - Y{i});
        end
    else
        % Use serial computation
        for i = 1:task_num
            grad_W(:, i) = X{i} * (X{i}' * W(:, i) - Y{i});
        end
    end
end

% smooth part gradient.
function [funcVal] = funVal_eval(W)
    funcVal = 0;
    
    % Use parallel computation if opts.pFlag is true
    if opts.pFlag
        parfor i = 1:task_num
            funcVal = funcVal + 0.5 * norm(Y{i} - X{i}' * W(:, i))^2;
        end
    else
        for i = 1:task_num
            funcVal = funcVal + 0.5 * norm(Y{i} - X{i}' * W(:, i))^2;
        end
    end
end

function [non_smooth_value] = nonsmooth_eval(W, rho_1, rho_2, rho_3)
    non_smooth_value = 0;
    
    % Pre-compute the size of W for efficiency
    num_rows = size(W, 1);
    
    for i = 1:num_rows
        w = W(i, :);
        non_smooth_value = non_smooth_value + ...
            rho_1 * norm(w, 1) + ...
            rho_2 * norm(R * w', 1) + ...
            rho_3 * norm(w, 2);
    end
end


end