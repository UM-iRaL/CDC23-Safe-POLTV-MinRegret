function [Phi_xw, Phi_uw, Phi_xe, Phi_ue, objective] = causal_constrained(sys, sls, opt, flag)
%CAUSAL_CONSTRAINED computes a constrained causal linear control policy
%that is optimal either in the H2 or in the Hinf sense

    % Define the decision variables of the optimization problem
    Phi_xw = sdpvar(sys.n*opt.T, sys.n*opt.T, 'full');
    Phi_uw = sdpvar(sys.m*opt.T, sys.n*opt.T, 'full');
    Phi_xe = sdpvar(sys.n*opt.T, sys.p*opt.T, 'full');
    Phi_ue = sdpvar(sys.m*opt.T, sys.p*opt.T, 'full');
    Phi_   = [Phi_xw Phi_xe; Phi_uw Phi_ue];
    
    % Define the objective function
    if strcmp(flag, 'H2')
        objective = norm(sqrtm(opt.C)*Phi_, 'fro');
    elseif strcmp(flag, 'Hinf')
        objective = norm(sqrtm(opt.C)*Phi_, 2);
    else
        error('Something went wrong...');
    end
    
    constraints = [];
    % Impose the achievability constraints
    constraints = [constraints, (sls.I - sls.Z*sls.A)*Phi_xw - sls.Z*sls.B*Phi_uw == sls.I];
    constraints = [constraints, (sls.I - sls.Z*sls.A)*Phi_xe - sls.Z*sls.B*Phi_ue == sls.Oxe];
    constraints = [constraints, Phi_xw*(sls.I - sls.Z*sls.A) - Phi_xe*sls.C == sls.I];
    constraints = [constraints, Phi_uw*(sls.I - sls.Z*sls.A) - Phi_ue*sls.C == sls.Ouw];

    % Impose the causal sparsities on the closed loop responses
    for i = 0:opt.T-2
        for j = i+1:opt.T-1 % Set j from i+2 for non-strictly causal controller (first element in w is x0)
            constraints = [constraints, Phi_xw((1+i*sys.n):((i+1)*sys.n), (1+j*sys.n):((j+1)*sys.n)) == zeros(sys.n, sys.n)];
            constraints = [constraints, Phi_uw((1+i*sys.m):((i+1)*sys.m), (1+j*sys.n):((j+1)*sys.n)) == zeros(sys.m, sys.n)];
            constraints = [constraints, Phi_xe((1+i*sys.n):((i+1)*sys.n), (1+j*sys.p):((j+1)*sys.p)) == zeros(sys.n, sys.p)];
            constraints = [constraints, Phi_ue((1+i*sys.m):((i+1)*sys.m), (1+j*sys.p):((j+1)*sys.p)) == zeros(sys.m, sys.p)];
        end
    end

    % Impose the polytopic safety constraints
    z = sdpvar(size(sls.Hnoise, 1), size(sls.H, 1), 'full'); % Define the dual variables
    for i = 1:size(sls.H, 1)
        constraints = [constraints, z(:, i)' * sls.hnoise <= sls.h(i)];
        constraints = [constraints, z(:, i) >= 0];
    end
    constraints = [constraints, sls.H * Phi_ == z' * sls.Hnoise];
    
    % Solve the optimization problem
    options = sdpsettings('verbose', 0, 'solver', 'mosek');
    sol = optimize(constraints, objective, options);
    if ~(sol.problem == 0)
        error('Something went wrong...');
    end
    
    % Extract the closed-loop responses corresponding to a constrained causal 
    % linear controller that is optimal either in the H2 or in the Hinf sense
    Phi_xw = value(Phi_xw); 
    Phi_uw = value(Phi_uw);
    Phi_xe = value(Phi_xe); 
    Phi_ue = value(Phi_ue);
    
    objective = value(objective)^2; % Extract the H2- or Hinf-optimal cost incurred by a constrained causal linear controller

end