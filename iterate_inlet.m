function [M_out, converged, gas_out] = iterate_inlet(mdot, A, gas_in, eta_i, M_in, gas_out)
    tol = 0.01; max_iter = 100; converged = false; n_iter = 0;
    
    a_in = sqrt(gas_in.gamma * gas_in.R * gas_in.T);
    V_in = M_in * a_in;
    T_0in = get_T(gas_in.T, gas_in.gamma, M_in);
    
    T_0out = T_0in;
    rho_in = gas_in.P / (gas_in.R * gas_in.T);
    V_out_guess = mdot / (rho_in * A);
    
    gas_out.gamma = gas_in.gamma;
    gas_out.cp = gas_in.cp;
    gas_out.R = gas_in.R;
    
    while ~converged && n_iter <= max_iter
        T_out = gas_in.T + (V_in^2 / (2 * gas_in.cp) - V_out_guess^2 / (2 * gas_out.cp));
        p_0out = gas_in.P * (1 + eta_i * V_in^2 / (2 * gas_in.cp * gas_in.T))^(gas_in.gamma / (gas_in.gamma - 1));
        p_out = p_0out * (T_out / T_0out)^(gas_out.gamma / (gas_out.gamma - 1));
        
        gas_out.T = T_out;
        gas_out.P = p_out;
        rho_out = p_out / (gas_out.R * T_out);
        
        V_out = V_out_guess;
        V_out_guess = mdot / (rho_out * A);
        
        if abs(V_out - V_out_guess) < tol
            converged = true;
            a_out = sqrt(gas_out.gamma * gas_out.R * gas_out.T);
            M_out = V_out / a_out;
        else
            n_iter = n_iter + 1;
        end
    end
    if ~converged, M_out = 0; end
end