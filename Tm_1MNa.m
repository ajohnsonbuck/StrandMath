function Tm = Tm_1MNa(dH0,dS0,conc)
    % Gas constant, cal/(mol K)
    R = 1.987204258; 
    
    % Calculate Tm from entropy, enthalpy, and concentration
    Tm = dH0/(dS0 + R*log(conc/4))-273.15;
end