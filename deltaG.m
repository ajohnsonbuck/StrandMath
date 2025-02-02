function dG = deltaG(dH0,dS0,T)
    T = T+273.15;
    dG = dH0 - T*dS0;
end