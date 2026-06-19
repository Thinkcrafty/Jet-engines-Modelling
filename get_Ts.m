function Ts = get_Ts(T0, gamma, M)
    Ts = T0 / (1 + ((gamma - 1)/2) * M^2);
end