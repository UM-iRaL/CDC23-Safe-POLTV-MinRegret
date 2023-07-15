clc; close all; clear;
addpath('./functions') % Add path to the folder with auxiliary functions
rng(1234);             % Set random seed for reproducibility

sys_num = 4;
control_num = 4;
horizon_max = 30;
noise_num = 8;

for sys_idx = 1:sys_num

    cost_mean = zeros(horizon_max-1, control_num, noise_num);  % store mean cost
    cost_std = zeros(horizon_max-1, control_num, noise_num);   % store standard deviation of cost

    %% Definition of the underlying discrete-time LTI system
    if sys_idx == 1 
        % % Open-loop stable system
        sys.rho = 0.85; % Spectral radius
        sys.A = sys.rho*[0.7 0.2 0; 0.3 0.7 -0.1; 0 -0.2 0.8];
        sys.B = [1 0.2; 2 0.3; 1.5 0.5];
        sys.C1 = [1 0 0; 0 1 0];
        sys.C2 = [0 1 0; 0 0 1];
        
    elseif sys_idx == 2
        % % Open-loop unstable system
        sys.rho = 1.05; % Spectral radius
        sys.A = sys.rho*[0.7 0.2 0; 0.3 0.7 -0.1; 0 -0.2 0.8];
        sys.B = [1 0.2; 2 0.3; 1.5 0.5];
        sys.C1 = [1 0 0; 0 1 0];
        sys.C2 = [0 1 0; 0 0 1];

    elseif sys_idx == 3
        sys.A = [1 0.15; 0 1];
        sys.B = [0.5; 0.5];
        sys.C1 = eye(2);
        sys.C2 = [1 0.5; 0.5 -1];

    else
        horizon_max = 25;
        [sys.A, sys.B] = UAV_model(9.81, 1);
        sys.C1 = [1 0 0 0 0 0;
                  0 1 0 0 0 0;
                  0 0 1 0 0 0];
        sys.C2 = [0 0 0 1 0 0;
                  0 0 0 0 1 0;
                  0 0 0 0 0 1];
    end

    sys.n = size(sys.A, 1);   % Order of the system: state dimension
    sys.m = size(sys.B, 2);   % Number of input channels
    sys.p = size(sys.C1, 1);   % Number of measurements
    sys.x0 = zeros(sys.n, 1); % Initial condition
    sys.noise_norm = 1;

    if sys_idx == 1
        sys.Hu = [eye(sys.m); -eye(sys.m)]; % Polytopic constraints: Hu * u <= hu
        sys.hu = 5*ones(size(sys.Hu, 1), 1);
        
        sys.Hx = [eye(sys.n); -eye(sys.n)]; % Polytopic constraints: Hx * x <= hx
        sys.hx = 5*ones(size(sys.Hx, 1), 1);
    end

    if sys_idx == 2 
        sys.Hu = [eye(sys.m); -eye(sys.m)]; % Polytopic constraints: Hu * u <= hu
        sys.hu = 30*ones(size(sys.Hu, 1), 1);
        
        sys.Hx = [eye(sys.n); -eye(sys.n)]; % Polytopic constraints: Hx * x <= hx
        sys.hx = 30*ones(size(sys.Hx, 1), 1);
    end

    if sys_idx == 3
        sys.Hu = [eye(sys.m); -eye(sys.m)]; % Polytopic constraints: Hu * u <= hu
        sys.hu = 30*ones(size(sys.Hu, 1), 1);
        
        sys.Hx = [eye(sys.n); -eye(sys.n)]; % Polytopic constraints: Hx * x <= hx
        sys.hx = 30*ones(size(sys.Hx, 1), 1);
    end

    if sys_idx == 4
        sys.Hu = [eye(sys.m); -eye(sys.m)]; % Polytopic constraints: Hu * u <= hu
        sys.hu = [pi; pi; 20; pi; pi; 20]; % 20*ones(size(sys.Hu, 1), 1);
        
        sys.Hx = [eye(sys.n); -eye(sys.n)]; % Polytopic constraints: Hx * x <= hx
        sys.hx = 5*ones(size(sys.Hx, 1), 1);   

        sys.noise_norm = 0.1;
    end
    
    sys.Hw = [eye(sys.n); -eye(sys.n)]; % Polytopic disturbance set: Hw * w <= hw 
    sys.hw = sys.noise_norm*ones(size(sys.Hw, 1), 1);
    
    sys.He = [eye(sys.p); -eye(sys.p)]; % Polytopic disturbance set: He * e <= he
    sys.he = sys.noise_norm*ones(size(sys.He, 1), 1);
    
    
    for T = 2:horizon_max

        h=waitbar(0, 'please wait');
        str = [num2str((T-1)/horizon_max*100) '%'];
        waitbar((T-1)/horizon_max, h, str)
        clear opt sls;

        %% Definition of the parameters of the optimization problem
        opt.Qt = eye(sys.n); % Stage cost: state weight matrix
        opt.Rt = eye(sys.m); % Stage cost: input weight matrix
        
        opt.T = T; % Control horizon
        
        opt.Q = kron(eye(opt.T), opt.Qt); % State cost matrix
        opt.R = kron(eye(opt.T), opt.Rt); % Input cost matrix
        opt.C = blkdiag(opt.Q, opt.R); % Cost matrix

        %% Definition of the stacked system dynamics over the control horizon
        sls.A = kron(eye(opt.T), sys.A);
        sls.B = kron(eye(opt.T), sys.B);
        if sys_idx < 4
            if mod(T,2) == 0
                sys.C = blkdiag(sys.C1, sys.C2);
                sls.C = kron(eye(opt.T/2), sys.C);
            else
                sys.C = blkdiag(sys.C1, sys.C2);
                sls.C = kron(eye(floor(opt.T/2)), sys.C);
                sls.C = blkdiag(sls.C, sys.C1);
            end
        else % for UAV, 1 position measurement followed by 2 velocity measurements
            if mod(T,3) == 0
                sys.C = blkdiag(sys.C1, sys.C2, sys.C2);
                sls.C = kron(eye(opt.T/3), sys.C);
            elseif mod(T,3) == 1
                sys.C = blkdiag(sys.C1, sys.C2, sys.C2);
                sls.C = kron(eye(floor(opt.T/3)), sys.C);
                sls.C = blkdiag(sls.C, sys.C1);
            else
                sys.C = blkdiag(sys.C1, sys.C2, sys.C2);
                sls.C = kron(eye(floor(opt.T/3)), sys.C);
                sls.C = blkdiag(sls.C, sys.C1, sys.C2);
            end    
        end

        sls.I = eye(sys.n*opt.T); % Identity matrix and block-downshift operator
        sls.Oxe = zeros(sys.n*opt.T, sys.p*opt.T);
        sls.Ouw = zeros(sys.m*opt.T, sys.n*opt.T);
        sls.Z = [zeros(sys.n, sys.n*(opt.T-1)) zeros(sys.n, sys.n); eye(sys.n*(opt.T-1)) zeros(sys.n*(opt.T-1), sys.n)];
        
        % Polytopic disturbance description and safety constraints
        sls.Hu = kron(eye(opt.T), sys.Hu);
        sls.hu = kron(ones(opt.T, 1), sys.hu);
        
        sls.Hx = kron(eye(opt.T), sys.Hx);
        sls.hx = kron(ones(opt.T, 1), sys.hx);
        
        sls.H = blkdiag(sls.Hu, sls.Hx);
        sls.h = [sls.hu; sls.hx];
        
        sls.Hw = kron(eye(opt.T), sys.Hw);
        sls.hw = kron(ones(opt.T, 1), sys.hw);
        
        sls.He = kron(eye(opt.T), sys.He);
        sls.he = kron(ones(opt.T, 1), sys.he);
        
        sls.Hnoise = blkdiag(sls.Hw, sls.He);
        sls.hnoise = [sls.hw; sls.he];

        %% Computation of the optimal noncausal unconstrained controller
        % The optimal dynamic sequence of control actions is unique.
        % However, the optimal H2 and Hinf costs incurred by the clairvoyant
        % controller are different!
        [Phi_nc_unc.xw, Phi_nc_unc.uw, Phi_nc_unc.xe, Phi_nc_unc.ue, obj_nc.unc_h2, obj_nc.unc_hinf] = noncausal_unconstrained(sys, sls, opt);
        % %% Computation of the optimal noncausal constrained H2 and Hinf controller
        % [Phi_nc_con_h2.xw, Phi_nc_con_h2.uw, Phi_nc_con_h2.xe, Phi_nc_con_h2.ue, obj_nc.con_h2]   = noncausal_constrained(sys, sls, opt, 'H2');
        % [Phi_nc_con_hinf.xw, Phi_nc_con_hinf.uw, Phi_nc_con_hinf.xe, Phi_nc_con_hinf.ue, obj_nc.con_hinf] = noncausal_constrained(sys, sls, opt, 'Hinf');
        % %% Computation of the optimal causal unconstrained H2 and Hinf controller
        % [Phi_c_unc_h2.xw, Phi_c_unc_h2.uw, Phi_c_unc_h2.xe, Phi_c_unc_h2.ue, obj_c.unc_h2]   = causal_unconstrained(sys, sls, opt, 'H2');
        % [Phi_c_unc_hinf.xw, Phi_c_unc_hinf.uw, Phi_c_unc_hinf.xe, Phi_c_unc_hinf.ue, obj_c.unc_hinf] = causal_unconstrained(sys, sls, opt, 'Hinf');
        %% Computation of the optimal causal constrained H2 and Hinf controller
        [Phi_c_con_h2.xw, Phi_c_con_h2.uw, Phi_c_con_h2.xe, Phi_c_con_h2.ue, obj_c.con_h2]   = causal_constrained(sys, sls, opt, 'H2');
        [Phi_c_con_hinf.xw, Phi_c_con_hinf.uw, Phi_c_con_hinf.xe, Phi_c_con_hinf.ue, obj_c.con_hinf] = causal_constrained(sys, sls, opt, 'Hinf');
        %% Computation of the regret-optimal causal controller 
        % [Phi_reg_unc_nc_unc.xw, Phi_reg_unc_nc_unc.uw, Phi_reg_unc_nc_unc.xe, Phi_reg_unc_nc_unc.ue, obj_reg.unc_nc_unc] = regret_unconstrained(sys, sls, opt, Phi_nc_unc);
        [Phi_reg_con_nc_unc.xw, Phi_reg_con_nc_unc.uw, Phi_reg_con_nc_unc.xe, Phi_reg_con_nc_unc.ue, obj_reg.con_nc_unc] = regret_constrained(sys, sls, opt, Phi_nc_unc);
        % [Phi_reg_con_nc_con_h2.xw, Phi_reg_con_nc_con_h2.uw, Phi_reg_con_nc_con_h2.xe, Phi_reg_con_nc_con_h2.ue, obj_reg.con_nc_con_h2]   = regret_constrained(sys, sls, opt, Phi_nc_con_h2);
        % [Phi_reg_con_nc_con_hinf.xw, Phi_reg_con_nc_con_hinf.uw, Phi_reg_con_nc_con_hinf.xe, Phi_reg_con_nc_con_hinf.ue, obj_reg.con_nc_con_hinf] = regret_constrained(sys, sls, opt, Phi_nc_con_hinf);

        %% Numerical experiments: comparison between time-averaged incurred control cost
        disturbance.profiles = ["Gaussian: N(0,1)" "Uniform: U(0.5, 1)" "Gamma: (4, 9)" "Exponential: mean 10" "Bernoulli: 0.5" "Weibull: (4, 9)" "Poisson: 4" "Worst-case"];
        disturbance.stochastic = [1000*ones(7, 1); ones(1, 1)]; % Maximum number of iterations per profile
        
        for i = 1:size(disturbance.profiles, 2) % Iterate over all different disturbance profiles
            for j = 1:disturbance.stochastic(i) 
                % Sample a disturbance realization
                if i == 1     % Gaussian: N(0, 1)
                    w_e_ = [sys.x0; randn(sys.n*(opt.T - 1), 1); randn(sys.p*opt.T, 1)];
                elseif i == 2 % Uniform: U(0, 1)
                    w_e_ = [sys.x0; rand(sys.n*(opt.T - 1), 1); rand(sys.p*opt.T, 1)];
                elseif i == 3 % Gamma: shape 4 and scale 9
                    w_e_ = [sys.x0; gamrnd(4,9,sys.n*(opt.T - 1), 1); gamrnd(4,9,sys.p*opt.T, 1)];
                elseif i == 4 % Exponential: mean 10 
                    w_e_ = [sys.x0; exprnd(10,sys.n*(opt.T - 1), 1); exprnd(10,sys.p*opt.T, 1)];
                elseif i == 5 % Bernoulli: 0.5
                    w_e_ = [sys.x0; binornd(1,0.5,sys.n*(opt.T - 1), 1); binornd(1,0.5,sys.p*opt.T, 1)];
                elseif i == 6 % Weibull: shape 4 and scale 9
                    w_e_ = [sys.x0; wblrnd(4,9,sys.n*(opt.T - 1), 1); wblrnd(4,9,sys.p*opt.T, 1)];
                elseif i == 7 % Poisson: rate parameter 4
                    w_e_ = [sys.x0; poissrnd(4,sys.n*(opt.T - 1), 1); poissrnd(4,sys.p*opt.T, 1)];
                else % Worst-case disturbance: adversarial selection for all three safe control laws 
                    % Compute the matrix that defines the quadratic form for the
                    % cost incurred by the H2, Hinfinity and regret-optimal controller
                    Phi_c_con_hinf_all = [Phi_c_con_hinf.xw, Phi_c_con_hinf.xe; Phi_c_con_hinf.uw, Phi_c_con_hinf.ue];
                    c_con_hinf.cost_qf = Phi_c_con_hinf_all' * opt.C * Phi_c_con_hinf_all;
                    
                    Phi_c_con_h2_all = [Phi_c_con_h2.xw, Phi_c_con_h2.xe; Phi_c_con_h2.uw, Phi_c_con_h2.ue];
                    c_con_h2.cost_qf   = Phi_c_con_h2_all' * opt.C * Phi_c_con_h2_all;
                    
                    Phi_reg_con_nc_unc_all = [Phi_reg_con_nc_unc.xw, Phi_reg_con_nc_unc.xe; Phi_reg_con_nc_unc.uw, Phi_reg_con_nc_unc.ue];
                    reg_con_nc_unc.cost_qf = Phi_reg_con_nc_unc_all' * opt.C * Phi_reg_con_nc_unc_all;

                    Phi_nc_unc_all = [Phi_nc_unc.xw, Phi_nc_unc.xe; Phi_nc_unc.uw, Phi_nc_unc.ue];
                    nc_unc.cost_qf = Phi_nc_unc_all' * opt.C * Phi_nc_unc_all;
                   
                    % Extract the direction of the eigenvector associated with the
                    % largest eigenvalue while maintaining zero initial condition
                    [c_con_hinf.evectors, c_con_hinf.evalues] = eig(c_con_hinf.cost_qf(sys.n+1:end, sys.n+1:end), 'vector');
                    [c_con_h2.evectors,   c_con_h2.evalues]   = eig(c_con_h2.cost_qf(sys.n+1:end, sys.n+1:end), 'vector');
                    [reg_con_nc_unc.evectors, reg_con_nc_unc.evalues] = eig(reg_con_nc_unc.cost_qf(sys.n+1:end, sys.n+1:end), 'vector');
                    [nc_unc.evectors, nc_unc.evalues] = eig(nc_unc.cost_qf(sys.n+1:end, sys.n+1:end), 'vector');
                    
                    [~, c_con_hinf.index] = max(c_con_hinf.evalues);
                    [~, c_con_h2.index]   = max(c_con_h2.evalues);
                    [~, reg_con_nc_unc.index]  = max(reg_con_nc_unc.evalues);
                    [~, nc_unc.index]   = max(nc_unc.evalues);

                    c_con_hinf.w_e_ = [sys.x0; c_con_hinf.evectors(:, c_con_hinf.index)]; 
                    c_con_h2.w_e_   = [sys.x0; c_con_h2.evectors(:, c_con_h2.index)]; 
                    reg_con_nc_unc.w_e_  = [sys.x0; reg_con_nc_unc.evectors(:, reg_con_nc_unc.index)]; 
                    nc_unc.w_e_  = [sys.x0; nc_unc.evectors(:, nc_unc.index)]; 
                end
                if i ~= 8 % Always simulate the three considered safe control policies with the same disturbance sequence, except when dealing with the worst-case of each of them
                    c_con_hinf.w_e_ = sys.noise_norm*w_e_/norm(w_e_); 
                    c_con_h2.w_e_ = sys.noise_norm*w_e_/norm(w_e_); 
                    reg_con_nc_unc.w_e_ = sys.noise_norm*w_e_/norm(w_e_);
                    nc_unc.w_e_ = sys.noise_norm*w_e_/norm(w_e_); 
                end
                % Vectorize the sampled disturbance sequence
                c_con_hinf.w_e_ = c_con_hinf.w_e_(:); 
                c_con_h2.w_e_ = c_con_h2.w_e_(:); 
                reg_con_nc_unc.w_e_ = reg_con_nc_unc.w_e_(:); 
                nc_unc.w_e_ = nc_unc.w_e_(:); 
                
                % Simulate the closed-loop system with the optimal causal constrained H2 and Hinf controller
                c_con_h2.cum_costs(j)   = evaluate_policy(opt, Phi_c_con_h2, c_con_h2.w_e_); 
                c_con_hinf.cum_costs(j) = evaluate_policy(opt, Phi_c_con_hinf, c_con_hinf.w_e_);
                % Simulate the closed-loop system with the regret-optimal causal controller
                reg_con_nc_unc.cum_costs(j) = evaluate_policy(opt, Phi_reg_con_nc_unc, reg_con_nc_unc.w_e_);
                % Simulate the closed-loop system with the noncausal constrained H2 controller
                nc_unc.cum_costs(j) = evaluate_policy(opt, Phi_nc_unc, nc_unc.w_e_);
            end
            
            cost_mean(T-1, :, i) = [mean(c_con_h2.cum_costs,'omitnan'), mean(c_con_hinf.cum_costs,'omitnan'), ...
                                    mean(reg_con_nc_unc.cum_costs,'omitnan'), mean(nc_unc.cum_costs,'omitnan')];
            cost_std(T-1, :, i) = [std(c_con_h2.cum_costs,'omitnan'), std(c_con_hinf.cum_costs,'omitnan'), ...
                                    std(reg_con_nc_unc.cum_costs,'omitnan'), std(nc_unc.cum_costs,'omitnan')];

            clear c_con_h2 c_con_hinf reg_con_nc_unc nc_unc; % Clear variables corresponding to past disturbances profiles
            
        end
        clear i j w_e_;
        delete(h);

    end
    filename = ['data_system_' num2str(sys_idx) '_T_1_' num2str(horizon_max)];
    save(filename, "cost_mean", "cost_std");

end