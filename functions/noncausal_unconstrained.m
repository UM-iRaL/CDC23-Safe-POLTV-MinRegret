function [Phi_xw, Phi_uw, Phi_xe, Phi_ue, obj_h2, obj_hinf] = noncausal_unconstrained(sys, sls, opt)
%NONCAUSAL_UNCONSTRAINED computes the unconstrained clairvoyant controller 
%that selects the globally optimal dynamic sequence of control actions in 
%hindsight
    
    % Define the decision variables of the optimization problem
    Phi_xw = sdpvar(sys.n*opt.T, sys.n*opt.T, 'full');
    Phi_uw = sdpvar(sys.m*opt.T, sys.n*opt.T, 'full');
%     Phi_xe = sdpvar(sys.n*opt.T, sys.p*opt.T, 'full');
%     Phi_ue = sdpvar(sys.m*opt.T, sys.p*opt.T, 'full');
    Phi_xe = zeros(sys.n*opt.T, sys.p*opt.T);
    Phi_ue = zeros(sys.m*opt.T, sys.p*opt.T);
    Phi_   = [Phi_xw Phi_xe; Phi_uw Phi_ue];

    % Define the objective function
    objective = norm(sqrtm(opt.C)*Phi_, 'fro');
    
    constraints = [];
    % Impose the achievability constraints
    constraints = [constraints, (sls.I - sls.Z*sls.A)*Phi_xw - sls.Z*sls.B*Phi_uw == sls.I];
%     constraints = [constraints, (sls.I - sls.Z*sls.A)*Phi_xe - sls.Z*sls.B*Phi_ue == sls.Oxe];
%     constraints = [constraints, Phi_xw*(sls.I - sls.Z*sls.A) - Phi_xe*sls.C == sls.I];
%     constraints = [constraints, Phi_uw*(sls.I - sls.Z*sls.A) - Phi_ue*sls.C == sls.Ouw];

    % Solve the optimization problem
    options = sdpsettings('verbose', 0, 'solver', 'mosek');
    sol = optimize(constraints, objective, options);
    if ~(sol.problem == 0)
        error('Something went wrong...');
    end
    
    % Extract the closed-loop responses corresponding to the unconstrained clairvoyant controller
    Phi_xw = value(Phi_xw); 
    Phi_uw = value(Phi_uw);
%     Phi_xe = value(Phi_xe); 
%     Phi_ue = value(Phi_ue);
    Phi_   = [Phi_xw Phi_xe; Phi_uw Phi_ue];

    obj_h2 = value(objective)^2;                       % Extract the H2-optimal   cost incurred by the unconstrained clairvoyant controller
    obj_hinf = norm(sqrtm(opt.C)*Phi_, 2)^2; % Compute the Hinf-optimal cost incurred by the unconstrained clairvoyant controller  
end