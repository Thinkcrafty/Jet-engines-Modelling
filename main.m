clear; clc; close all;

%% 02 - Engine Parameters Definitions
eng_param = struct();
eng_perf  = struct();

% Inlet
eng_param.A1 = 0.30;       % [m^2] inlet area at station 1
eng_param.A2 = 0.32;       % [m^2] inlet exit area at comp face
eng_perf.eta_i = 0.98;     % inlet efficiency

% Compressor
eng_param.comp_n_stages = 10;
eng_perf.CPR = 6.1;        % overall compressor pressure rise
eng_perf.eta_c = 0.80;     % isentropic stage efficiency

% Combustor
eng_perf.eta_b = 0.90;     % combustion efficiency
eng_perf.dp_over_p = 0.06; % total pressure loss (%)
eng_perf.max_f = 0.25;     % max fuel (percent of Stoichiometric)
eng_perf.min_f = 0.125;    % min fuel 
eng_perf.V_nominal = 45;   % [m/s] nominal flow velocity
eng_perf.LHV = 43.2e6;     % [J/kg] Jet fuel lower heating value (approx)

% Turbine
eng_param.turb_n_stages = 2;
eng_perf.eta_t = 0.80;     % isentropic stage efficiency
eng_perf.mech_loss = 0.99; % mechanical loss

% Nozzle
eng_param.A8 = 0.27;       % [m^2] nozzle area
eng_perf.eta_noz = 0.80;   % nozzle losses

%% 04 - Engine Ambient / Operating Conditions
eng_op_con = struct();
eng_op_con.throttle_pos = 1.0; 
eng_op_con.mdot_guess = 20;     % [kg/s]
eng_op_con.alt = 35000;         % [ft]
eng_op_con.M_i = 0.8;

M_i = eng_op_con.M_i;
V_i = M2Vt(M_i, eng_op_con.alt) * 0.514444; % True airspeed [m/s] (kts to m/s)
p_amb = isa_p(eng_op_con.alt);  % Static pressure [Pa]
T_amb = isa_T(eng_op_con.alt);  % Static temperature [K]

%% 05 - Initial Stations Setup
station_names = {'ambient', 'inlet', 'inlet @ comp face', 'after compressor', 'after combustor', 'after turbine', 'nozzle exit'};

gas = repmat(struct('T', 0, 'P', 0, 'gamma', 1.4, 'cp', 1004, 'R', 287.05), 1, 7);
M = zeros(1, 7);

% Initialize ambient station
for i = 1:7
    gas(i).T = T_amb;
    gas(i).P = p_amb;
    M(i) = M_i;
end

%% 06 - From Ambient (a) to Station 1
[M_calc, conv, gas(2)] = iterate_inlet(eng_op_con.mdot_guess, eng_param.A1, gas(1), 1.0, M(1), gas(2));
if conv, M(2) = M_calc; end

%% 07 - From Station 1 to Station 2
[M_calc, conv, gas(3)] = iterate_inlet(eng_op_con.mdot_guess, eng_param.A2, gas(2), eng_perf.eta_i, M(2), gas(3));
if conv
    M(3:end) = M_calc; % Assume constant Mach through machine initially
else
    disp('ERROR: inlet did not converge');
end

%% 09 - Compressor Exit Station (2 -> 3)
[stages_T_out, stages_p_out, conv, comp_work, gas(4)] = multi_stage_compressor(gas(3), eng_param.comp_n_stages, eng_perf.CPR, eng_perf.eta_c, M(3), gas(4));

%% 11 - Combustor (3 -> 4)

phi = (eng_perf.max_f - eng_perf.min_f) * eng_op_con.throttle_pos + eng_perf.min_f;
f_stoich = 0.068; % approx Stoichiometric ratio for kerosene
f = phi * f_stoich;
mixt_frac = f / (1 + f);

[M_calc, conv, gas(5)] = iterate_combustor(gas(4), eng_perf.V_nominal, M(4), eng_perf.dp_over_p, gas(5), f, eng_perf.LHV, eng_perf.eta_b);
if conv
    M(5) = M_calc;
else
    disp('ERROR: combustor did not converge');
end

% Propagate gas composition downstream (hot section)
for i = 6:7
    gas(i).gamma = 1.33;
    gas(i).cp = 1156;
end

%% 13 - Turbine (4 -> 5)
[turb_T_out, turb_p_out, turb_work, gas(6)] = multi_stage_turbine(gas(5), comp_work, eng_param.turb_n_stages, eng_perf.eta_t, eng_perf.mech_loss, M(5), M(6), gas(6));

%% 15 - Nozzle (5 -> 8)
[choked, mdot_noz, M(7), F, gas(7)] = calc_nozzle(gas(6), M(6), eng_perf.eta_noz, p_amb, eng_param.A8, V_i, gas(7));

%% Output Results
fprintf('=== TURBOJET CYCLE ANALYSIS ===\n');
for i = 1:7
    fprintf('Station %d (%-18s): Mach %5.3f | p = %5.2f atm, T = %5.0f K, p0 = %5.2f atm, T0 = %5.0f K\n', ...
        i, station_names{i}, M(i), gas(i).P/101325, gas(i).T, ...
        get_p(gas(i).P, gas(i).gamma, M(i))/101325, get_T(gas(i).T, gas(i).gamma, M(i)));
end

fprintf('\n=== PERFORMANCE METRICS ===\n');
fprintf('Delta mdot iter : %.2f kg/s\n', mdot_noz - eng_op_con.mdot_guess);
fprintf('Thrust          : %.0f kN\n', (F * mdot_noz)/1000);
fprintf('Alt             : %.0f ft (T_amb = %.1f K, Mach = %.2f)\n', eng_op_con.alt, T_amb, M_i);
fprintf('Mass Flow       : %.0f kg/s\n', mdot_noz);

mdot_fuel = (mixt_frac / eng_perf.eta_b) * mdot_noz;
TSFC = (mdot_fuel / mdot_noz) / F;
fprintf('Fuel Flow       : %.0f kg/h\n', mdot_fuel * 3600);
fprintf('TSFC            : %.2f kg/(kN h)\n', TSFC * 3600 * 1000);

SAR = V_i / mdot_fuel;
fprintf('Spec. Air Range : %.1f m/kg (%.4f nm/kg)\n', SAR, V_i * 1.94384 / (mdot_fuel * 3600));
