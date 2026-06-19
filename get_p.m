function p0 = get_p(ps, gamma, M)
    p0 = ps * (1 + ((gamma - 1)/2) * M^2)^(gamma / (gamma - 1));
end
