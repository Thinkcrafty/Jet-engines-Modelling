function ps = get_ps(p0, Ts, T0, gamma)
    ps = p0 * (Ts / T0)^(gamma / (gamma - 1));
end