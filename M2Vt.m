function Vt = M2Vt(M, Hc_ft)
    gamma = 1.4; R = 287.053; ms2kt = 1 / 0.514444;
    Vt = M * sqrt(gamma * R * isa_T(Hc_ft)) * ms2kt; % Returns kts
end