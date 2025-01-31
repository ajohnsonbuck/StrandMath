% Salt correction of oligonucleotide hybridization as a function of salt
% (Na+ and Mg2+) concentration.
%
% Based on salt corrections used on IDT OligoAnalyzer, from Owczarzy,R.
% et al., Biochemistry, 43, 3537 and Owczarzy, R. et al., Biochemistry, 47, 5336
% Tm_1MNa = melting temperature (deg C) at 1 M Na+
% fGC = fractional GC content
% Na = [Na+], in molar (moles/L)

function Tm = salt_correction(Tm_1MNa, Nbp, fGC, Na, Mg)

    % Convert from Celsius to Kelvin
    Tm_1MNa = Tm_1MNa + 273.15;
    
    R = sqrt(Mg)/Na;
    if R < 0.22
        Tm = 1/(1/Tm_1MNa + ((4.29*fGC - 3.95)*log(Na) + 0.940*(log(Na))^2)*(1E-5)); % Monovalent correction
    elseif R < 6
        a = 3.92*(0.843-0.352*sqrt(Na)*log(Na));
        d = 1.42*(1.279-4.03E-3*log(Na)-8.03E-3*(log(Na))^2);
        g = 8.31*(0.486-0.258*log(Na)+5.25E-3*(log(Na))^3);
        Tm = 1/(1/Tm_1MNa + (a - 0.911*log(Mg)+fGC*(6.26+d*log(Mg))+1/(2*(Nbp-1))*(-48.2+52.5*log(Mg)+g*(log(Mg))^2))*1E-5);
    else
        a = 3.92;
        d = 1.42;
        g = 8.31;
        Tm = 1/(1/Tm_1MNa + (a - 0.911*log(Mg)+fGC*(6.26+d*log(Mg)+1/(2*(Nbp-1))*(-48.2+52.5*log(Mg)+g*log(Mg))))*1E-5);
    end

    % Convert from Kelvin back to Celsius
    Tm = Tm-273.15;

end