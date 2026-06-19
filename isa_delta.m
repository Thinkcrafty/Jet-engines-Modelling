function d = isa_delta(Hc_ft)
    L = -0.0065; T_SL = 288.15; R = 287.053; g_SL = 9.80665;
    if Hc_ft <= 36089.24
        d = (1 + (L / T_SL) * (Hc_ft / 3.28084))^(-g_SL / (L * R));
    else
        d = 0.22336 * exp(-(g_SL / (R * 216.65)) * ((Hc_ft - 36089.25) / 3.28084));
    end
end