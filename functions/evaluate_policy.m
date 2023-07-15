function [cum_cost] = evaluate_policy(opt, Phi, noise)
%EVALUATE_POLICY computes the cumulative cost incurred applying the policy 
%corresponding to the closed-loop responses in Phi in response to 
%the disturbance realization w, e
    
    % Compute the input-state trajectory associated with the disturbance w 
    Phi_all = [Phi.xw Phi.xe; Phi.uw Phi.ue];
    x_u_ = Phi_all * noise;

    % Compute the incurred cumulative cost
    cum_cost = x_u_'*opt.C*x_u_;
    
end