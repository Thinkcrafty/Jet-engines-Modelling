function [stages_T, stages_p, converged, comp_work, gas_out] = multi_stage_compressor(gas_in, n_stages, CPR, eta_c, M_in, gas_out)
T_0in_orig = get_T(gas_in.T, gas_in.gamma, M_in);
p_0in_orig = get_p(gas_in.P, gas_in.gamma, M_in);

CR_stage = CPR^(1 / n_stages);
stage_mult = ones(1, n_stages);
shifter = ones(1, n_stages);
shift = 0.001;
step_shift = shift / n_stages;

center = floor(n_stages / 2);
for i = 1:center
    shifter(i) = 1 + (center - i + 1) * step_shift;
    shifter(n_stages - i + 1) = 1 - ((center - i + 1) * step_shift);
end

converged = false; n_iter = 0; max_iter = 5000; prev_delta_t = 1000;
stages_p = zeros(1, n_stages);
stages_T = zeros(1, n_stages);

gas_out = gas_in; % Init properties

while ~converged && n_iter <= max_iter
    T_0in = T_0in_orig;
    p_0in = p_0in_orig;
    comp_work = 0;
    stage_gas = gas_in;

    for st = 1:n_stages
        T_i = stage_gas.T;
        p0 = p_0in * CR_stage * stage_mult(st);
        T0 = T_0in / eta_c * ((p0 / p_0in)^((stage_gas.gamma - 1) / stage_gas.gamma) - 1) + T_0in;

        T = get_Ts(T0, stage_gas.gamma, M_in);
        p = get_ps(p0, T, T0, stage_gas.gamma);

        stage_gas.T = T;
        stage_gas.P = p;

        stages_p(st) = p;
        stages_T(st) = T;
        comp_work = comp_work + stage_gas.cp * (T - T_i);

        p_0in = p0;
        T_0in = T0;
    end

    if n_stages > 2
        max_delta_t = max(diff(stages_T));
    elseif n_stages > 1
        max_delta_t = max([stages_T(2) - stages_T(1), stages_T(1) - gas_in.T]);
    else
        max_delta_t = T - gas_in.T;
    end

    if max_delta_t < prev_delta_t && n_iter < max_iter
        n_iter = n_iter + 1;
        stage_mult = stage_mult .* shifter;
        prev_delta_t = max_delta_t;
    elseif n_iter >= max_iter
        n_iter = n_iter + 1;
    else
        converged = true;
    end
end
gas_out.T = T;
gas_out.P = p;
end