function Tm = Tm_1MNa(dH0,dS0,conc)
    R = 1.987204258; % cal/(mol K)

    Tm = dH0/(dS0 + R*log(conc))-273.15;
end