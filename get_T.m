function T0 = get_T(Ts, gamma, M)
    T0 = Ts * (1 + ((gamma - 1)/2) * M^2);
end