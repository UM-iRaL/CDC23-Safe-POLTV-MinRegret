function [x, u] = evaluate_state_input(opt, Phi, noise)
%EVALUATE_STATE_INPUT computes the cost of state and input incurred applying the policy 
%corresponding to the closed-loop responses in Phi in response to 
%the disturbance realization w, e
    
    Phi_x = [Phi.xw Phi.xe]; 
    x = opt.Q * (Phi_x * noise) .* (Phi_x * noise);

    Phi_u = [Phi.uw Phi.ue];
    u = opt.R * (Phi_u * noise) .* (Phi_u * noise);
    
end