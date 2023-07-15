function [x, u] = evaluate_traj(opt, sys, Phi, noise)
%EVALUATE_STATE_INPUT computes the state and input incurred applying the policy 
%corresponding to the closed-loop responses in Phi in response to 
%the disturbance realization w, e
    
    Phi_x = [Phi.xw Phi.xe]; 
    x = Phi_x * noise;
    x = reshape(x, sys.n, opt.T);

    Phi_u = [Phi.uw Phi.ue];
    u = Phi_u * noise;
    u = reshape(u, sys.m, opt.T);
end