function [M_out, converged, gas_out] = iterate_combustor(gas_in, V_nom, M_in, dp_p, gas_out, f, LHV, eta_b)
tol = 0.01; max_iter = 100; converged = false; n_iter = 0;

T_0in = get_T(gas_in.T, gas_in.gamma, M_in);
p_0in = get_p(gas_in.P, gas_in.gamma, M_in);

% Thermodynamic balance to replace gas.equilibrate('HP')
% (1) * cp_cold * T_0in + f * LHV * eta_b = (1+f) * cp_hot * T_0out
gas_out.gamma = 1.33; % hot gas prop
gas_out.cp = 1156;    % hot gas prop
gas_out.R = 287.05;

T_0out = (gas_in.cp * T_0in + f * LHV * eta_b) / ((1 + f) * gas_out.cp);
p_0out = p_0in * (1 - dp_p);

T_out = T_0out; % Initial guess
p_out = p_0out;

while ~converged && n_iter <= max_iter
    gas_out.T = T_out;
    gas_out.P = p_out;

    a_out = sqrt(gas_out.gamma * gas_out.R * gas_out.T);
    M_out = V_nom / a_out;

    T_out_new = get_Ts(T_0out, gas_out.gamma, M_out);
    p_out_new = get_ps(p_0out, T_out_new, T_0out, gas_out.gamma);

    if abs(gas_out.P - p_out_new) < tol
        converged = true;
        gas_out.T = T_out_new;
        gas_out.P = p_out_new;
    else
        T_out = T_out_new;
        p_out = p_out_new;
        n_iter = n_iter + 1;
    end
end
if ~converged, M_out = 0; end
end