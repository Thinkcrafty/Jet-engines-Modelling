function t = isa_theta(Hc_ft)
    if Hc_ft <= 36089.24
        t = (1 + (-0.0065 / 288.15) * (Hc_ft / 3.28084));
    else
        t = 216.65 / 288.15;
    end
end