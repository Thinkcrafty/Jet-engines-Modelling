function [stages_T, stages_p, turb_work, gas_out] = multi_stage_turbine(gas_in, W_c, n_stages, eta_t, eta_m, M_in, M_out, gas_out)
T_0in = get_T(gas_in.T, gas_in.gamma, M_in);
p_0in = get_p(gas_in.P, gas_in.gamma, M_in);

T_in = get_Ts(T_0in, gas_in.gamma, M_out);
p_in = get_ps(p_0in, T_in, T_0in, gas_in.gamma);

W_per_stage = (W_c / eta_m) / n_stages;
stages_p = zeros(1, n_stages);
stages_T = zeros(1, n_stages);

stage_gas = gas_in;
stage_gas.T = T_in;
stage_gas.P = p_in;
turb_work = 0;

gas_out = stage_gas; % Copy specific heats

for st = 1:n_stages
    T_i = stage_gas.T;
    T_0out_prime = T_0in - W_per_stage / (stage_gas.cp * eta_t);
    p_0out = p_0in * (T_0out_prime / T_0in)^(stage_gas.gamma / (stage_gas.gamma - 1));
    T_0out = T_0in - eta_t * (T_0in - T_0out_prime);

    T = get_Ts(T_0out, stage_gas.gamma, M_out);
    p = get_ps(p_0out, T, T_0out, stage_gas.gamma);

    stage_gas.T = T;
    stage_gas.P = p;

    stages_p(st) = p;
    stages_T(st) = T;
    turb_work = turb_work + stage_gas.cp * (T - T_i);

    p_0in = p_0out;
    T_0in = T_0out;
end
gas_out.T = T;
gas_out.P = p;
end
