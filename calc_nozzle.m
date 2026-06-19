function [choked, mdot, M_out, F, gas_out] = calc_nozzle(gas_in, M_in, eta_noz, p_amb, A_star, V_i, gas_out)
p0_in = get_p(gas_in.P, gas_in.gamma, M_in);
T0_in = get_T(gas_in.T, gas_in.gamma, M_in);

pc_ratio = 1 / (1 - (1 / eta_noz) * ((gas_in.gamma - 1)/(gas_in.gamma + 1)))^(gas_in.gamma / (gas_in.gamma - 1));
gas_out = gas_in;

if p_amb <= p0_in / pc_ratio
    choked = true;
    p = p0_in / pc_ratio;
    T = T0_in / ((gas_in.gamma + 1) / 2);
    V = sqrt(gas_in.gamma * gas_in.R * T);
    rho = p / (gas_in.R * T);
    mdot = rho * V * A_star;
    F = (V - V_i) + (A_star / mdot) * (p - p_amb);
    M_out = 1;
else
    choked = false;
    p = p_amb;
    T = T0_in - eta_noz * T0_in * (1 - (1 / (p0_in / p_amb)^((gas_in.gamma - 1) / gas_in.gamma)));
    V = sqrt(2 * gas_in.cp * (T0_in - T));
    rho = p / (gas_in.R * T);
    mdot = rho * V * A_star;
    F = (V - V_i);
    a = sqrt(gas_in.gamma * gas_in.R * T);
    M_out = V / a;
end
gas_out.T = T;
gas_out.P = p;
end