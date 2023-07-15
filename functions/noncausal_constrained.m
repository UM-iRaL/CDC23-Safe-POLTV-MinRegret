function [Phi_xw, Phi_uw, Phi_xe, Phi_ue, objective] = noncausal_constrained(sys, sls, opt, flag)
%NONCAUSAL_CONSTRAINED computes a constrained clairvoyant linear control 
%policy that is optimal either in the H2 or in the Hinf sense
    
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
    
    % Extract the closed-loop responses corresponding to a constrained
    % clairvoyant linear controller that is optimal either in the H2 or in the Hinf sense
    Phi_xw = value(Phi_xw); 
    Phi_uw = value(Phi_uw);
    Phi_xe = value(Phi_xe); 
    Phi_ue = value(Phi_ue);
    
    objective = value(objective)^2; % Extract the H2- or Hinf-optimal cost incurred by a constrained clairvoyant linear controller

end